/*
Copyright (c) 2014-2015 Xiaowei Zhu, Tsinghua University

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <string>
#include <thread>
#include <vector>
#include <limits>

#include "core/atomic.hpp"
#include "core/constants.hpp"
#include "core/filesystem.hpp"
#include "core/partition.hpp"
#include "core/queue.hpp"
#include "core/time.hpp"
#include "core/type.hpp"

#include "GrReader.h"

long PAGESIZE = 4096;

void generate_edge_grid(std::string input, std::string output,
                        VertexId vertices, int partitions, int edge_type) {
  int parallelism = std::thread::hardware_concurrency();
  int edge_unit;

  ////////////////////////////////////////////////////////////////////////////////
  // TODO look here
  ////////////////////////////////////////////////////////////////////////////////
  // update edge unit AND edges
  EdgeId edges;
  switch (edge_type) {
  case 0:
    edge_unit = sizeof(VertexId) * 2;
    break;
  case 1:
    printf("FILLING EVERYTHING WITH DEFAULT WEIGHT OF 1\n");
    edge_unit = sizeof(VertexId) * 2 + sizeof(Weight);
    break;
  default:
    fprintf(stderr, "edge type (%d) is not supported.\n", edge_type);
    exit(-1);
  }
  printf("vertices = %d\n", vertices);
  ////////////////////////////////////////////////////////////////////////////////
  // TODO end here
  ////////////////////////////////////////////////////////////////////////////////

  char **buffers = new char *[parallelism * 2];
  bool *occupied = new bool[parallelism * 2];
  for (int i = 0; i < parallelism * 2; i++) {
    buffers[i] = (char *)memalign(PAGESIZE, IOSIZE);
    occupied[i] = false;
  }
  Queue<std::tuple<int, long>> tasks(parallelism);
  int **fout;
  std::mutex **mutexes;
  fout = new int *[partitions];
  mutexes = new std::mutex *[partitions];
  if (file_exists(output)) {
    remove_directory(output);
  }
  create_directory(output);

  const int grid_buffer_size = 768; // 12 * 8 * 8
  char *global_grid_buffer =
      (char *)memalign(PAGESIZE, grid_buffer_size * partitions * partitions);
  char ***grid_buffer = new char **[partitions];
  int **grid_buffer_offset = new int *[partitions];
  for (int i = 0; i < partitions; i++) {
    mutexes[i] = new std::mutex[partitions];
    fout[i] = new int[partitions];
    grid_buffer[i] = new char *[partitions];
    grid_buffer_offset[i] = new int[partitions];
    for (int j = 0; j < partitions; j++) {
      char filename[4096];
      sprintf(filename, "%s/block-%d-%d", output.c_str(), i, j);
      fout[i][j] = open(filename, O_WRONLY | O_APPEND | O_CREAT, 0644);
      grid_buffer[i][j] =
          global_grid_buffer + (i * partitions + j) * grid_buffer_size;
      grid_buffer_offset[i][j] = 0;
    }
  }

  std::vector<std::thread> threads;
  for (int ti = 0; ti < parallelism; ti++) {
    threads.emplace_back([&]() {
      char *local_buffer = (char *)memalign(PAGESIZE, IOSIZE);
      int *local_grid_offset = new int[partitions * partitions];
      int *local_grid_cursor = new int[partitions * partitions];
      VertexId source, target;
      Weight weight;
      while (true) {
        int cursor;
        long bytes;
        std::tie(cursor, bytes) = tasks.pop();
        if (cursor == -1)
          break;
        memset(local_grid_offset, 0, sizeof(int) * partitions * partitions);
        memset(local_grid_cursor, 0, sizeof(int) * partitions * partitions);
        char *buffer = buffers[cursor];
        for (long pos = 0; pos < bytes; pos += edge_unit) {
          source = *(VertexId *)(buffer + pos);
          target = *(VertexId *)(buffer + pos + sizeof(VertexId));
          int i = get_partition_id(vertices, partitions, source);
          int j = get_partition_id(vertices, partitions, target);
          local_grid_offset[i * partitions + j] += edge_unit;
        }
        local_grid_cursor[0] = 0;
        for (int ij = 1; ij < partitions * partitions; ij++) {
          local_grid_cursor[ij] = local_grid_offset[ij - 1];
          local_grid_offset[ij] += local_grid_cursor[ij];
        }
        assert(local_grid_offset[partitions * partitions - 1] == bytes);
        for (long pos = 0; pos < bytes; pos += edge_unit) {
          source = *(VertexId *)(buffer + pos);
          target = *(VertexId *)(buffer + pos + sizeof(VertexId));
          int i = get_partition_id(vertices, partitions, source);
          int j = get_partition_id(vertices, partitions, target);
          *(VertexId *)(local_buffer + local_grid_cursor[i * partitions + j]) =
              source;
          *(VertexId *)(local_buffer + local_grid_cursor[i * partitions + j] +
                        sizeof(VertexId)) = target;
          if (edge_type == 1) {
            weight = *(Weight *)(buffer + pos + sizeof(VertexId) * 2);
            *(Weight *)(local_buffer + local_grid_cursor[i * partitions + j] +
                        sizeof(VertexId) * 2) = weight;
          }
          local_grid_cursor[i * partitions + j] += edge_unit;
        }
        int start = 0;
        for (int ij = 0; ij < partitions * partitions; ij++) {
          assert(local_grid_cursor[ij] == local_grid_offset[ij]);
          int i = ij / partitions;
          int j = ij % partitions;
          std::unique_lock<std::mutex> lock(mutexes[i][j]);
          if (local_grid_offset[ij] - start > edge_unit) {
            write(fout[i][j], local_buffer + start,
                  local_grid_offset[ij] - start);
          } else if (local_grid_offset[ij] - start == edge_unit) {
            memcpy(grid_buffer[i][j] + grid_buffer_offset[i][j],
                   local_buffer + start, edge_unit);
            grid_buffer_offset[i][j] += edge_unit;
            if (grid_buffer_offset[i][j] == grid_buffer_size) {
              write(fout[i][j], grid_buffer[i][j], grid_buffer_size);
              grid_buffer_offset[i][j] = 0;
            }
          }
          start = local_grid_offset[ij];
        }
        occupied[cursor] = false;
      }
    });
  }

  // open file descriptor
  int fin = -1;

  // TODO inaccurate measure of things
  //long total_bytes = file_size(input);
  // tracks progress (assuming all bytes are known)
  //long read_bytes = 0;

  // graph processing timer
  double start_time = get_time();

  // TODO currently will not read edge data
  galois::graphs::GrReader<void> graphReader;

  // split graphs into this amount of chunks to read (to prevent
  // loading entire graph at once into memory)
  uint64_t GRAPH_CHUNKS = 64;
  assert(GRAPH_CHUNKS > 0);
  uint64_t curChunk = 0;
  // read chunk by chunk
  graphReader.loadGraphBalanceEdges(input, curChunk, GRAPH_CHUNKS);

  uint64_t globalNodes = graphReader.size();
  uint64_t globalEdges = graphReader.sizeEdges();

  edges = globalEdges;

  // first, determine if this graph being read fits within GridGraph type
  // parameters
  uint64_t maxVertexID = std::numeric_limits<VertexId>::max();
  uint64_t maxEdgeID = std::numeric_limits<EdgeId>::max();
  //printf("Max vertex ID is %lu\n", maxVertexID);
  //printf("Max edge ID is %lu\n", maxEdgeID);
  if ((globalNodes > maxVertexID) || (globalEdges > maxEdgeID)) {
    fprintf(stderr, "Graph vertex/edge ID exceeds supported GG limit\n");
    exit(-1);
  }

  // total edges read
  uint64_t edgesRead = 0;

  // nodes/edges in this chunk to read
  uint64_t chunkEdgesRead = 0;
  uint64_t edgesInChunk = graphReader.initEdgeIterator();

  // determine how many edges to read per round via IOSIZE
  uint64_t maxEdgesToRead = IOSIZE / edge_unit;
  //printf("Reading %lu edges per chunk\n", maxEdgesToRead);

  // signifies which location a thread should read (cursor is pushed as a task)
  int cursor = 0;
  uint64_t totalBytes = 0;
  Weight dummyWeight = 1;

  // loop until all edges read
  while (edgesRead < globalEdges) {
    if (chunkEdgesRead == edgesInChunk) {
      // load next chunk if this chunk has no edges left
      curChunk++;
      graphReader.resetAndFree();
      graphReader.loadGraphBalanceEdges(input, curChunk, GRAPH_CHUNKS);
      chunkEdgesRead = 0;
      edgesInChunk = graphReader.initEdgeIterator();
      continue;
    }

    //////////////////////////////////////////////////////////////////////////
    // if code makes it this far, it means the current chunk I have loaded
    // has edges to read
    //////////////////////////////////////////////////////////////////////////

    // determine how many edges to read in this iteration; min of what's
    // left and IOSIZE
    uint64_t edgesLeftInChunk = edgesInChunk - chunkEdgesRead;
    uint64_t edgesToReadThisIteration = std::min(edgesLeftInChunk,
                                                 maxEdgesToRead);

    char* curBuffer = buffers[cursor];
    uint64_t bytesWritten = 0;

    // this loop will read the appropriate amount of edges
    for (uint64_t i = 0; i < edgesToReadThisIteration; i++) {
      uint64_t src;
      uint64_t dst;
      // read edge
      std::tie(src, dst) = graphReader.nextEdge();
      // "cast" to vertex id type
      VertexId cSrc = src;
      VertexId cDst = dst;
      //fprintf(stderr, "%lu %d %lu %d\n", src, cSrc, dst, cDst);

      // memcpy src and dst over to buffer
      memcpy(curBuffer + bytesWritten, &cSrc, sizeof(cSrc));
      memcpy(curBuffer + bytesWritten + sizeof(cSrc), &cDst, sizeof(cDst));
      bytesWritten += sizeof(VertexId) * 2;

      if (edge_type) {
        memcpy(curBuffer + bytesWritten, &dummyWeight, sizeof(Weight));
        bytesWritten += sizeof(Weight);
      }

      // below is to make sure it wrote to buffer correctly
      //VertexId* test = (VertexId*)curBuffer;
      //printf("%d %d\n", test[2 * i], test[2 * i + 1]);
    }

    totalBytes += bytesWritten;

    chunkEdgesRead += edgesToReadThisIteration;
    edgesRead += edgesToReadThisIteration;

    // start a thread on processing that read buffer
    occupied[cursor] = true;
    tasks.push(std::make_tuple(cursor, bytesWritten));
    printf("progress: %.2f%%\r", 100. * edgesRead / globalEdges);
    fflush(stdout);

    // find an unoccupied buffer to read the next thing into
    while (occupied[cursor]) {
      cursor = (cursor + 1) % (parallelism * 2);
    }
  }

  assert(edgesRead == globalEdges);

  // pushing kill signals
  for (int ti = 0; ti < parallelism; ti++) {
    tasks.push(std::make_tuple(-1, 0));
  }

  // start tasks here I guess
  for (int ti = 0; ti < parallelism; ti++) {
    threads[ti].join();
  }

  printf("%lf -> ", get_time() - start_time);
  long ts = 0;
  for (int i = 0; i < partitions; i++) {
    for (int j = 0; j < partitions; j++) {
      if (grid_buffer_offset[i][j] > 0) {
        ts += grid_buffer_offset[i][j];
        write(fout[i][j], grid_buffer[i][j], grid_buffer_offset[i][j]);
      }
    }
  }
  printf("%lf (%ld)\n", get_time() - start_time, ts);

  for (int i = 0; i < partitions; i++) {
    for (int j = 0; j < partitions; j++) {
      close(fout[i][j]);
    }
  }

  printf("it takes %.2f seconds to generate edge blocks\n",
         get_time() - start_time);

  long offset;
  int fout_column =
      open((output + "/column").c_str(), O_WRONLY | O_APPEND | O_CREAT, 0644);
  int fout_column_offset = open((output + "/column_offset").c_str(),
                                O_WRONLY | O_APPEND | O_CREAT, 0644);
  offset = 0;
  // writing
  for (int j = 0; j < partitions; j++) {
    for (int i = 0; i < partitions; i++) {
      printf("progress: %.2f%%\r", 100. * offset / totalBytes);
      fflush(stdout);
      write(fout_column_offset, &offset, sizeof(offset));
      char filename[4096];
      sprintf(filename, "%s/block-%d-%d", output.c_str(), i, j);
      offset += file_size(filename);
      fin = open(filename, O_RDONLY);
      while (true) {
        long bytes = read(fin, buffers[0], IOSIZE);
        assert(bytes != -1);
        if (bytes == 0)
          break;
        write(fout_column, buffers[0], bytes);
      }
      close(fin);
    }
  }

  write(fout_column_offset, &offset, sizeof(offset));
  close(fout_column_offset);
  close(fout_column);
  printf("column oriented grid generated\n");
  int fout_row =
      open((output + "/row").c_str(), O_WRONLY | O_APPEND | O_CREAT, 0644);
  int fout_row_offset = open((output + "/row_offset").c_str(),
                             O_WRONLY | O_APPEND | O_CREAT, 0644);

  // writing
  offset = 0;
  for (int i = 0; i < partitions; i++) {
    for (int j = 0; j < partitions; j++) {
      printf("progress: %.2f%%\r", 100. * offset / totalBytes);
      fflush(stdout);
      write(fout_row_offset, &offset, sizeof(offset));
      char filename[4096];
      sprintf(filename, "%s/block-%d-%d", output.c_str(), i, j);
      offset += file_size(filename);
      fin = open(filename, O_RDONLY);
      while (true) {
        long bytes = read(fin, buffers[0], IOSIZE);
        assert(bytes != -1);
        if (bytes == 0)
          break;
        write(fout_row, buffers[0], bytes);
      }
      close(fin);
    }
  }

  write(fout_row_offset, &offset, sizeof(offset));
  close(fout_row_offset);
  close(fout_row);
  printf("row oriented grid generated\n");

  printf("it takes %.2f seconds to generate edge grid\n",
         get_time() - start_time);

  FILE *fmeta = fopen((output + "/meta").c_str(), "w");
  fprintf(fmeta, "%d %d %ld %d", edge_type, vertices, edges, partitions);
  fclose(fmeta);
}

int main(int argc, char **argv) {
  int opt;
  std::string input = "";
  std::string output = "";
  VertexId vertices = -1;
  int partitions = -1;
  int edge_type = 0;
  while ((opt = getopt(argc, argv, "i:o:v:p:t:")) != -1) {
    switch (opt) {
    case 'i':
      input = optarg;
      break;
    case 'o':
      output = optarg;
      break;
    case 'v':
      vertices = atoi(optarg);
      break;
    case 'p':
      partitions = atoi(optarg);
      break;
    case 't':
      edge_type = atoi(optarg);
      break;
    }
  }
  if (input == "" || output == "" || vertices == -1) {
    fprintf(stderr,
            "usage: %s -i [input path] -o [output path] -v [vertices] -p "
            "[partitions] -t [edge type: 0=unweighted, 1=weighted]\n",
            argv[0]);
    exit(-1);
  }
  if (partitions == -1) {
    partitions = vertices / CHUNKSIZE;
  }
  generate_edge_grid(input, output, vertices, partitions, edge_type);
  return 0;
}
