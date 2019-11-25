/**
 * @file GrReader.h
 *
 * Contains the implementation of GrReader
 */

#ifndef __GR_READER__
#define __GR_READER__

#include "stdlib.h"
#include "unistd.h"
#include <fstream>
#include <utility>
#include <iostream>

namespace galois {
namespace graphs {

/**
 * Class to get access to a Galois GRs prefix sum without actually reading it
 * into memory.
 */
class OfflinePrefixSum {
  std::ifstream fileEdgeDst;
  std::streamoff locEdgeDst;

  uint64_t numNodes;
  uint64_t numEdges;
  uint64_t sizeEdgeData;
  uint64_t length;
  bool v2;

  uint64_t outIndexs(uint64_t node) {
    std::streamoff pos = (4 + node) * sizeof(uint64_t);
    // move to correct position in file
    if (locEdgeDst != pos) {
      fileEdgeDst.seekg(pos, fileEdgeDst.beg);
      locEdgeDst = pos;
    }

    // read the value
    uint64_t retval;
    try {
      fileEdgeDst.read(reinterpret_cast<char*>(&retval), sizeof(uint64_t));
    } catch (std::ifstream::failure e) {
      std::cerr << "Exception while reading edge destinations:" << e.what()
                << "\n";
      std::cerr << "IO error flags: EOF " << fileEdgeDst.eof() << " FAIL "
                << fileEdgeDst.fail() << " BAD " << fileEdgeDst.bad() << "\n";
    }

    // metadata update
    auto numBytesRead = fileEdgeDst.gcount();
    assert(numBytesRead == sizeof(uint64_t));
    locEdgeDst += numBytesRead;

    return retval;
  }

public:
  typedef uint64_t GraphNode;

  OfflinePrefixSum(const std::string& name)
      : fileEdgeDst(name, std::ios_base::binary), locEdgeDst(0) {
    if (!fileEdgeDst.is_open() || !fileEdgeDst.good()) throw "Bad filename";

    fileEdgeDst.exceptions(std::ifstream::eofbit | std::ifstream::failbit |
                           std::ifstream::badbit);
    uint64_t ver = 0;

    try {
      fileEdgeDst.read(reinterpret_cast<char*>(&ver), sizeof(uint64_t));
      fileEdgeDst.read(reinterpret_cast<char*>(&sizeEdgeData),
                       sizeof(uint64_t));
      fileEdgeDst.read(reinterpret_cast<char*>(&numNodes), sizeof(uint64_t));
      fileEdgeDst.read(reinterpret_cast<char*>(&numEdges), sizeof(uint64_t));
    } catch (std::ifstream::failure e) {
      std::cerr << "Exception while reading graph header:" << e.what() << "\n";
      std::cerr << "IO error flags: EOF " << fileEdgeDst.eof() << " FAIL "
                << fileEdgeDst.fail() << " BAD " << fileEdgeDst.bad() << "\n";
    }

    if (ver == 0 || ver > 2) throw "Bad Version";
    v2 = ver == 2;
    if (!fileEdgeDst) throw "Out of data";

    // File length
    fileEdgeDst.seekg(0, fileEdgeDst.end);
    length = fileEdgeDst.tellg();
    if (length < sizeof(uint64_t) * (4 + numNodes) +
                     (v2 ? sizeof(uint64_t) : sizeof(uint32_t)) * numEdges)
      throw "File too small";

    fileEdgeDst.seekg(0, std::ios_base::beg);
  }

  OfflinePrefixSum(OfflinePrefixSum&&) = default;

  /**
   * Accesses the prefix sum on disk.
   *
   * @param n Index into edge prefix sum
   * @returns The value located at index n in the edge prefix sum array
   */
  uint64_t operator[](uint64_t n) { return outIndexs(n); }
};

/**
 * Class that loads a portion of a Galois graph from disk directly into
 * memory buffers for access.
 *
 * @tparam EdgeDataType type of the edge data
 * @todo version 2 Galois binary graph support; currently only suppports
 * version 1
 */
template <typename EdgeDataType>
class GrReader {
private:
  // buffers that you load data into
  //! buffer that tells you where a particular node's edges begin
  uint64_t* outIndexBuffer = nullptr;
  //! buffer that tells the destination of edges
  uint32_t* edgeDestBuffer = nullptr;
  //! buffer that stores edge data
  EdgeDataType* edgeDataBuffer = nullptr;

  //! size of the entire graph (not just locallly loaded portion)
  uint32_t globalSize = 0;
  //! number of edges in the entire graph (not just locallly loaded portion)
  uint64_t globalEdgeSize = 0;

  //! number of nodes loaded into this graph
  uint32_t numLocalNodes = 0;
  //! number of edges loaded into this graph
  uint64_t numLocalEdges = 0;

  //! specifies how many nodes are skipped from the beginning of the graph
  //! in this loaded portion of it
  uint64_t nodeOffset = 0;
  //! specifies how many edges are skipped from the beginning of the graph
  //! in this loaded portion of it
  uint64_t edgeOffset = 0;
  //! specifies whether or not the graph is loaded
  bool graphLoaded = false;
  //! bool specifying if read file had edge data to load
  bool noEdgeDataExists = true;

  // accumulators for tracking bytes read

  /**
   * Load the out indices (i.e. where a particular node's edges begin in the
   * array of edges) from the file.
   *
   * @param graphFile loaded file for the graph
   * @param nodeStart the first node to load
   * @param numNodesToLoad number of nodes to load
   */
  void loadOutIndex(std::ifstream& graphFile, uint64_t nodeStart,
                    uint64_t numNodesToLoad) {
    if (numNodesToLoad == 0) {
      return;
    }
    assert(outIndexBuffer == nullptr);
    outIndexBuffer = (uint64_t*)malloc(sizeof(uint64_t) * numNodesToLoad);

    if (outIndexBuffer == nullptr) {
      fprintf(stderr, "cannot allocate mem for out index\n");
      exit(1);
    }

    // position to start of contiguous chunk of nodes to read
    uint64_t readPosition = (4 + nodeStart) * sizeof(uint64_t);
    graphFile.seekg(readPosition);

    uint64_t numBytesToLoad = numNodesToLoad * sizeof(uint64_t);
    uint64_t bytesRead      = 0;

    while (numBytesToLoad > 0) {
      graphFile.read(((char*)this->outIndexBuffer) + bytesRead, numBytesToLoad);
      size_t numRead = graphFile.gcount();
      numBytesToLoad -= numRead;
      bytesRead += numRead;
    }

    assert(numBytesToLoad == 0);

    nodeOffset = nodeStart;
  }

  /**
   * Load the edge destination information from the file.
   *
   * @param graphFile loaded file for the graph
   * @param edgeStart the first edge to load
   * @param numEdgesToLoad number of edges to load
   * @param numGlobalNodes total number of nodes in the graph file; needed
   * to determine offset into the file
   */
  void loadEdgeDest(std::ifstream& graphFile, uint64_t edgeStart,
                    uint64_t numEdgesToLoad, uint64_t numGlobalNodes) {
    if (numEdgesToLoad == 0) {
      return;
    }

    assert(edgeDestBuffer == nullptr);
    edgeDestBuffer = (uint32_t*)malloc(sizeof(uint32_t) * numEdgesToLoad);

    if (edgeDestBuffer == nullptr) {
      fprintf(stderr, "cannot allocate mem for edge dest\n");
      exit(1);
    }

    // position to start of contiguous chunk of edges to read
    uint64_t readPosition = (4 + numGlobalNodes) * sizeof(uint64_t) +
                            (sizeof(uint32_t) * edgeStart);
    graphFile.seekg(readPosition);

    uint64_t numBytesToLoad = numEdgesToLoad * sizeof(uint32_t);
    uint64_t bytesRead      = 0;
    while (numBytesToLoad > 0) {
      graphFile.read(((char*)this->edgeDestBuffer) + bytesRead, numBytesToLoad);
      size_t numRead = graphFile.gcount();
      numBytesToLoad -= numRead;
      bytesRead += numRead;
    }

    assert(numBytesToLoad == 0);
    // save edge offset of this graph for later use
    edgeOffset = edgeStart;
  }

  /**
   * Load the edge data information from the file.
   *
   * @tparam EdgeType must be non-void in order to call this function
   *
   * @param edgeStart the first edge to load
   * @param numEdgesToLoad number of edges to load
   * @param numGlobalNodes total number of nodes in the graph file; needed
   * to determine offset into the file
   * @param numGlobalEdges total number of edges in the graph file; needed
   * to determine offset into the file
   */
  template <
      typename EdgeType,
      typename std::enable_if<!std::is_void<EdgeType>::value>::type* = nullptr>
  void loadEdgeData(std::ifstream& graphFile, uint64_t edgeStart,
                    uint64_t numEdgesToLoad, uint64_t numGlobalNodes,
                    uint64_t numGlobalEdges) {
    if (numEdgesToLoad == 0) {
      return;
    }

    // if file has nothing to read then read nothing and return
    if (noEdgeDataExists) {
      return;
    }

    assert(edgeDataBuffer == nullptr);
    edgeDataBuffer =
        (EdgeDataType*)malloc(sizeof(EdgeDataType) * numEdgesToLoad);

    if (edgeDataBuffer == nullptr) {
      fprintf(stderr, "cannot allocate mem for edge data\n");
      exit(1);
    }

    // position after nodes + edges
    uint64_t baseReadPosition = (4 + numGlobalNodes) * sizeof(uint64_t) +
                                (sizeof(uint32_t) * numGlobalEdges);

    // version 1 padding TODO make version agnostic
    if (numGlobalEdges % 2) {
      baseReadPosition += sizeof(uint32_t);
    }

    // jump to first byte of edge data
    uint64_t readPosition =
        baseReadPosition + (sizeof(EdgeDataType) * edgeStart);
    graphFile.seekg(readPosition);
    uint64_t numBytesToLoad = numEdgesToLoad * sizeof(EdgeDataType);
    uint64_t bytesRead      = 0;

    while (numBytesToLoad > 0) {
      graphFile.read(((char*)this->edgeDataBuffer) + bytesRead, numBytesToLoad);
      size_t numRead = graphFile.gcount();
      numBytesToLoad -= numRead;
      bytesRead += numRead;
    }

    assert(numBytesToLoad == 0);
  }

  /**
   * Load edge data function for when the edge data type is void, i.e.
   * no edge data to load.
   *
   * Does nothing of importance.
   *
   * @tparam EdgeType if EdgeType is void, this function will be used
   */
  template <
      typename EdgeType,
      typename std::enable_if<std::is_void<EdgeType>::value>::type* = nullptr>
  void loadEdgeData(std::ifstream& graphFile, uint64_t edgeStart,
                    uint64_t numEdgesToLoad, uint64_t numGlobalNodes,
                    uint64_t numGlobalEdges) {
    // do nothing (edge data is void, i.e. no edge data)
  }

  /**
   * Resets graph metadata to default values. Does NOT touch the buffers.
   */
  void resetGraphStatus() {
    graphLoaded    = false;
    globalSize     = 0;
    globalEdgeSize = 0;
    nodeOffset     = 0;
    edgeOffset     = 0;
    numLocalNodes  = 0;
    numLocalEdges  = 0;
  }

  /**
   * Free all of the buffers in memory.
   */
  void freeMemory() {
    free(outIndexBuffer);
    outIndexBuffer = nullptr;
    free(edgeDestBuffer);
    edgeDestBuffer = nullptr;
    free(edgeDataBuffer);
    edgeDataBuffer = nullptr;
  }

  // TODO documentation
  // TODO can be template like original block range
  std::pair<uint32_t, uint32_t> blockRange(uint32_t begin, uint32_t end,
                                           unsigned id, unsigned num) {
    uint32_t dist   = end - begin;
    uint32_t numper = std::max((dist + num - 1) / num, (uint32_t)1); // round up
    uint32_t A      = std::min(numper * id, dist);
    uint32_t B      = std::min(numper * (id + 1), dist);
    begin += A;
    if (dist != B) {
      end = begin;
      end += (B - A);
    }
    return std::make_pair(begin, end);
  }

  //! Helper function that fills up a vector with consecutive numbers
  uint32_t setupScaleFactor(uint32_t numDivisions,
                            std::vector<unsigned>& scaleFactor) {
    uint32_t numBlocks = 0;
    // if scale factor isn't specified, everyone gets the same amount
    numBlocks = numDivisions;

    // scale factor holds a prefix sum of the scale factor
    for (uint32_t i = 0; i < numDivisions; i++) {
      scaleFactor.push_back(i + 1);
    }
    return numBlocks;
  }

  //! Given a lower and upper bound, try to find a point that matches the
  //! target weight given a prefix sum.
  template <typename PrefixSumType>
  size_t findIndexPrefixSum(size_t targetWeight, uint64_t lb, uint64_t ub,
                            PrefixSumType& edgePrefixSum) {
    while (lb < ub) {
      size_t mid = lb + (ub - lb) / 2;
      size_t num_edges;

      if (mid != 0) {
        num_edges = edgePrefixSum[mid - 1];
      } else {
        num_edges = 0;
      }

      size_t weight = num_edges + mid;

      if (weight < targetWeight) {
        lb = mid + 1;
      } else if (weight >= targetWeight) {
        ub = mid;
      }
    }

    return lb;
  }

  //! Gives user a node split that attempts to balance edge counts.
  template <typename NodeType = uint64_t>
  std::pair<NodeType, NodeType> divideNodesBinarySearch(NodeType numNodes,
      uint64_t numEdges, size_t id, size_t total, std::string fileToRead) {
    //using NodeRange = std::pair<uint32_t, uint32_t>;

    // numNodes = 0 corner case
    if (numNodes == 0) {
      return std::make_pair((NodeType)0, (NodeType)0);
    }

    assert(total >= 1);
    assert(id >= 0 && id < total);

    // prep scale factor
    std::vector<unsigned> scaleFactor;
    scaleFactor.reserve(total);
    // scale factor setup
    uint32_t numBlocks = setupScaleFactor(total, scaleFactor);
    // weight of all data
    uint64_t weight = numEdges + 1;
    //uint64_t weight = numEdges;
    // weight of a block (one block for each division by default; if scale
    // factor specifies something different, then use that instead)
    uint64_t blockWeight = (weight + numBlocks - 1) / numBlocks;
    //printf("block weight %lu\n", blockWeight);

    // lower and upper blocks that this division should use determined
    // using scaleFactor
    uint32_t blockLower;
    if (id != 0) {
      blockLower = scaleFactor[id - 1];
    } else {
      blockLower = 0;
    }
    uint32_t blockUpper = scaleFactor[id];
    assert(blockLower <= blockUpper);

    // offline prefix sum class
    OfflinePrefixSum pSum(fileToRead);
    //printf("[%d] block lower is %lu, block upper is %lu\n", id, blockLower, blockUpper);

    // lower and upper
    uint64_t nodesLower;
    // use prefix sum to find node bounds
    if (blockLower == 0) {
      nodesLower = 0;
    } else {
      nodesLower = findIndexPrefixSum(blockWeight * blockLower, 0, numNodes,
                                      pSum);
    }
    uint64_t nodesUpper;
    nodesUpper = findIndexPrefixSum(blockWeight * blockUpper, nodesLower,
                                    numNodes, pSum);

    if (id == (total - 1)) {
      if (nodesUpper != numNodes) {
        printf("[%d] Readjusting last host to get all nodes\n", id);
        nodesUpper = numNodes;
      }
    }

    if (id == 0) {
      printf("Total nodes %lu, Total edges %lu\n", numNodes, numEdges);
    }
    printf("[%d] Lower is %lu, upper is %lu\n", id, nodesLower, nodesUpper);
    return std::make_pair(nodesLower, nodesUpper);
  }

public:
  /**
   * Class vars should be initialized by in-class initialization; all
   * left is to reset read counters.
   */
  GrReader() { }

  /**
   * On destruction, free allocated buffers (if necessary).
   */
  ~GrReader() noexcept { freeMemory(); }

  // copy not allowed
  //! disabled copy constructor
  GrReader(const GrReader&) = delete;
  //! disabled copy constructor operator
  GrReader& operator=(const GrReader&) = delete;
  // move not allowed
  //! disabled move operator
  GrReader(GrReader&&) = delete;
  //! disabled move constructor operator
  GrReader& operator=(GrReader&&) = delete;

  /**
   * Gets the number of global nodes in the graph
   * @returns the total number of nodes in the graph (not just local loaded
   * nodes)
   */
  uint32_t size() const { return globalSize; }

  /**
   * Gets the number of local nodes in the graph
   */
  uint32_t sizeLocal() const { return numLocalNodes; }

  /**
   * Gets the number of global edges in the graph
   * @returns the total number of edges in the graph (not just local loaded
   * edges)
   */
  uint64_t sizeEdges() const { return globalEdgeSize; }

  /**
   * Gets the number of local edges in the graph
   */
  uint32_t sizeEdgesLocal() const { return numLocalEdges; }

  //! @returns node offset of this buffered graph
  uint64_t getNodeOffset() const { return nodeOffset; }

  /**
   * Loads given Galois CSR graph into memory.
   *
   * @param filename name of graph to load; should be in Galois binary graph
   * format
   */
  void loadGraph(const std::string& filename) {
    if (graphLoaded) {
      fprintf(stderr, "cannot load buffered graph more than once\n");
      exit(1);
    }

    std::ifstream graphFile(filename.c_str());
    uint64_t header[4];
    graphFile.read(((char*)header), sizeof(uint64_t) * 4);

    numLocalNodes = globalSize = header[2];
    numLocalEdges = globalEdgeSize = header[3];

    loadOutIndex(graphFile, 0, globalSize);
    loadEdgeDest(graphFile, 0, globalEdgeSize, globalSize);
    // may or may not do something depending on EdgeDataType
    loadEdgeData<EdgeDataType>(graphFile, 0, globalEdgeSize, globalSize,
                               globalEdgeSize);
    graphLoaded = true;

    graphFile.close();
  }

  /**
   * Read in a chunk of a graph based on given ID; balance edges read
   * by each host.
   *
   * Some hosts may read nothing.
   *
   * TODO finish this documentation
   * TODO reduce code duplication
   */
  void loadGraphBalanceEdges(const std::string& filename, uint32_t myID,
                             uint32_t totalIDs) {
    // TODO
    if (graphLoaded) {
      fprintf(stderr, "cannot load graph more than once\n");
      exit(1);
    }

    std::ifstream graphFile(filename.c_str());
    graphFile.seekg(0);
    // read header
    uint64_t header[4];
    graphFile.read(((char*)header), sizeof(uint64_t) * 4);

    if (header[1] == 0) {
      printf("No edge data exists to read from graph! Default to 0 weights.\n");
      noEdgeDataExists = true;
    } else {
      assert(header[1] == 4);
      noEdgeDataExists = false;
    }

    globalSize     = header[2];
    globalEdgeSize = header[3];

    // split nodes among all participants
    std::pair<uint32_t, uint32_t> myRange =
      divideNodesBinarySearch(globalSize, globalEdgeSize, myID, totalIDs,
                              filename);
    uint32_t nodeStart = myRange.first;
    uint32_t nodeEnd = myRange.second;
    //printf("start is %lu, end is %lu\n", nodeStart, nodeEnd);

    // first load out indices
    numLocalNodes = nodeEnd - nodeStart;
    loadOutIndex(graphFile, nodeStart, numLocalNodes);

    // get the edge offset
    OfflinePrefixSum pSum(filename);
    if (nodeStart == 0) {
      edgeOffset = 0;
    } else {
      edgeOffset = pSum[nodeStart - 1];
    }

    // set this early so i can use edgebegin/end
    graphLoaded = true;

    uint64_t edgeStart;
    uint64_t lastEdge;

    // determine how many edges we will read using loaded out indices
    if (numLocalNodes > 0) {
      edgeStart = edgeBegin(nodeStart);
      lastEdge = edgeEnd(nodeEnd - 1);
    } else {
      numLocalEdges = 0;
      edgeStart = 0;
      lastEdge = 0;
    }
    numLocalEdges = lastEdge - edgeStart;

    //printf("edge start is %lu, last edge is %lu\n", edgeStart, lastEdge);

    // load destinations and data
    loadEdgeDest(graphFile, edgeStart, numLocalEdges, globalSize);
    loadEdgeData<EdgeDataType>(graphFile, edgeStart, numLocalEdges,
                               globalSize, globalEdgeSize);
    // finish up
    graphFile.close();
  }


  /**
   * Read in a chunk of a graph based on given ID; balance nodes read
   * by each host.
   *
   * Some hosts may read nothing.
   *
   * TODO finish this documentation
   */
  void loadGraphBalanceNodes(const std::string& filename, uint32_t myID,
                             uint32_t totalIDs) {
    // TODO
    if (graphLoaded) {
      fprintf(stderr, "cannot load graph more than once\n");
      exit(1);
    }

    std::ifstream graphFile(filename.c_str());
    graphFile.seekg(0);
    // read header
    uint64_t header[4];
    graphFile.read(((char*)header), sizeof(uint64_t) * 4);

    if (header[1] == 0) {
      printf("No edge data exists to read from graph! Default to 0 weights.");
      noEdgeDataExists = true;
    } else {
      assert(header[1] == 4);
      noEdgeDataExists = false;
    }

    globalSize     = header[2];
    globalEdgeSize = header[3];

    // split nodes among all participants
    std::pair<uint32_t, uint32_t> myRange = blockRange(0, globalSize, myID,
                                                       totalIDs);
    uint32_t nodeStart = myRange.first;
    uint32_t nodeEnd = myRange.second;

    // first load out indices
    numLocalNodes = nodeEnd - nodeStart;
    loadOutIndex(graphFile, nodeStart, numLocalNodes);

    // get the edge offset
    OfflinePrefixSum pSum(filename);
    if (nodeStart == 0) {
      edgeOffset = 0;
    } else {
      edgeOffset = pSum[nodeStart - 1];
    }

    // set this early so i can use edgebegin/end
    graphLoaded = true;

    uint64_t edgeStart;
    uint64_t lastEdge;

    // determine how many edges we will read using loaded out indices
    if (numLocalNodes > 0) {
      edgeStart = edgeBegin(nodeStart);
      lastEdge = edgeEnd(nodeEnd - 1);
    } else {
      numLocalEdges = 0;
      edgeStart = 0;
      lastEdge = 0;
    }
    numLocalEdges = lastEdge - edgeStart;


    // load destinations and data
    loadEdgeDest(graphFile, edgeStart, numLocalEdges, globalSize);
    loadEdgeData<EdgeDataType>(graphFile, edgeStart, numLocalEdges,
                               globalSize, globalEdgeSize);
    // finish up
    graphFile.close();
  }


  /**
   * Given a node/edge range to load, loads the specified portion of the graph
   * into memory buffers using read.
   *
   * @param filename name of graph to load; should be in Galois binary graph
   * format
   * @param nodeStart First node to load
   * @param nodeEnd Last node to load, non-inclusive
   * @param edgeStart First edge to load; should correspond to first edge of
   * first node
   * @param edgeEnd Last edge to load, non-inclusive
   * @param numGlobalNodes Total number of nodes in the graph
   * @param numGlobalEdges Total number of edges in the graph
   */
  void loadPartialGraph(const std::string& filename, uint64_t nodeStart,
                        uint64_t nodeEnd, uint64_t edgeStart, uint64_t edgeEnd,
                        uint64_t numGlobalNodes, uint64_t numGlobalEdges) {
    if (graphLoaded) {
      fprintf(stderr, "cannot load graph more than once\n");
      exit(1);
    }

    std::ifstream graphFile(filename.c_str());
    graphFile.seekg(0);
    // read header
    uint64_t header[4];
    graphFile.read(((char*)header), sizeof(uint64_t) * 4);

    if (header[1] == 0) {
      printf("No edge data exists to read from graph! Default to 0 weights.");
      noEdgeDataExists = true;
    } else {
      assert(header[1] == 4);
      noEdgeDataExists = false;
    }

    globalSize     = numGlobalNodes;
    globalEdgeSize = numGlobalEdges;

    assert(nodeEnd >= nodeStart);
    numLocalNodes = nodeEnd - nodeStart;
    loadOutIndex(graphFile, nodeStart, numLocalNodes);

    assert(edgeEnd >= edgeStart);
    numLocalEdges = edgeEnd - edgeStart;
    loadEdgeDest(graphFile, edgeStart, numLocalEdges, numGlobalNodes);

    // may or may not do something depending on EdgeDataType
    loadEdgeData<EdgeDataType>(graphFile, edgeStart, numLocalEdges,
                               numGlobalNodes, numGlobalEdges);
    graphLoaded = true;
    graphFile.close();
  }

  /**
   * Get the index to the first edge of the provided node THAT THIS GRAPH
   * HAS LOADED (not necessary the first edge of it globally).
   *
   * @param globalNodeID the global node id of the node to get the edge
   * for
   * @returns a GLOBAL edge id
   */
  uint64_t edgeBegin(uint64_t globalNodeID) {
    if (!graphLoaded) {
      fprintf(stderr, "graph not loaded yet\n");
      exit(1);
    }

    if (numLocalNodes == 0) {
      return 0;
    }
    assert(nodeOffset <= globalNodeID);
    assert(globalNodeID < (nodeOffset + numLocalNodes));

    uint64_t localNodeID = globalNodeID - nodeOffset;

    if (localNodeID != 0) {
      return outIndexBuffer[localNodeID - 1];
    } else {
      return edgeOffset;
    }
  }

  /**
   * Get the index to the first edge of the node after the provided node.
   *
   * @param globalNodeID the global node id of the node to get the edge
   * for
   * @returns a GLOBAL edge id iterator
   */
  uint64_t edgeEnd(uint64_t globalNodeID) {
    if (!graphLoaded) {
      fprintf(stderr, "Graph not loaded yet\n");
      exit(1);
    }

    if (numLocalNodes == 0) {
      return 0;
    }
    assert(nodeOffset <= globalNodeID);
    assert(globalNodeID < (nodeOffset + numLocalNodes));

    uint64_t localNodeID = globalNodeID - nodeOffset;
    return outIndexBuffer[localNodeID];
  }

  /**
   * Get the global node id of the destination of the provided edge.
   *
   * @param globalEdgeID the global edge id of the edge to get the destination
   * for (should obtain from edgeBegin/End)
   */
  uint64_t edgeDestination(uint64_t globalEdgeID) {
    if (!graphLoaded) {
      fprintf(stderr, "Graph not loaded yet\n");
      exit(1);
    }

    if (numLocalEdges == 0) {
      return 0;
    }
    assert(edgeOffset <= globalEdgeID);
    assert(globalEdgeID < (edgeOffset + numLocalEdges));

    uint64_t localEdgeID = globalEdgeID - edgeOffset;
    return edgeDestBuffer[localEdgeID];
  }

  /**
   * Get the edge data of some edge.
   *
   * @param globalEdgeID the global edge id of the edge to get the data of
   * @returns the edge data of the requested edge id
   */
  template <typename K = EdgeDataType,
            typename std::enable_if<!std::is_void<K>::value>::type* = nullptr>
  EdgeDataType edgeData(uint64_t globalEdgeID) {
    if (!graphLoaded) {
      fprintf(stderr, "Graph not loaded yet\n");
      exit(1);
    }

    // return 0 for edge data if none exists
    if (noEdgeDataExists) {
      return 0;
    }

    // if exists but buffer is null, die
    if (edgeDataBuffer == nullptr) {
      fprintf(stderr, "getting edge data when graph has no edge data\n");
      exit(1);
    }

    if (numLocalEdges == 0) {
      return 0;
    }

    assert(edgeOffset <= globalEdgeID);
    assert(globalEdgeID < (edgeOffset + numLocalEdges));

    uint64_t localEdgeID = globalEdgeID - edgeOffset;
    return edgeDataBuffer[localEdgeID];
  }

  /**
   * Version of above function when edge data type is void.
   */
  template <typename K = EdgeDataType,
            typename std::enable_if<std::is_void<K>::value>::type* = nullptr>
  unsigned edgeData(uint64_t globalEdgeID) {
    return 0;
  }

  /**
   * Free all of the in memory buffers in this object and reset graph status.
   */
  void resetAndFree() {
    freeMemory();
    resetGraphStatus();
  }
};
} // namespace graphs
} // namespace galois
#endif
