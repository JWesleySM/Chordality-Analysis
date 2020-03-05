//===-- RegAllocBasic.cpp - Basic Register Allocator ----------------------===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
//
// This file defines the RABasic function pass, which provides a minimal
// implementation of the basic register allocator.
//
//===----------------------------------------------------------------------===//

#include "AllocationOrder.h"
#include "LiveDebugVariables.h"
#include "RegAllocBase.h"
#include "Spiller.h"
#include "llvm/Analysis/AliasAnalysis.h"
#include "llvm/CodeGen/CalcSpillWeights.h"
#include "llvm/CodeGen/LiveIntervals.h"
#include "llvm/CodeGen/LiveRangeEdit.h"
#include "llvm/CodeGen/LiveRegMatrix.h"
#include "llvm/CodeGen/LiveStacks.h"
#include "llvm/CodeGen/MachineBlockFrequencyInfo.h"
#include "llvm/CodeGen/MachineFunctionPass.h"
#include "llvm/CodeGen/MachineInstr.h"
#include "llvm/CodeGen/MachineLoopInfo.h"
#include "llvm/CodeGen/MachineRegisterInfo.h"
#include "llvm/CodeGen/Passes.h"
#include "llvm/CodeGen/RegAllocRegistry.h"
#include "llvm/CodeGen/TargetRegisterInfo.h"
#include "llvm/CodeGen/VirtRegMap.h"
#include "llvm/PassAnalysisSupport.h"
#include "llvm/Support/Debug.h"
#include "llvm/Support/raw_ostream.h"
#include <cstdlib>
#include <queue>

//JW includes
//Library necesary to get information about the Pysical register
#include "llvm/Target/TargetMachine.h"

using namespace llvm;

#define DEBUG_TYPE "regalloc"

static RegisterRegAlloc basicRegAlloc("basic", "basic register allocator",
                                      createBasicRegisterAllocator);

namespace {
  struct CompSpillWeight {
    bool operator()(LiveInterval *A, LiveInterval *B) const {
      return A->weight < B->weight;
    }
  };
}

namespace {
/// RABasic provides a minimal implementation of the basic register allocation
/// algorithm. It prioritizes live virtual registers by spill weight and spills
/// whenever a register is unavailable. This is not practical in production but
/// provides a useful baseline both for measuring other allocators and comparing
/// the speed of the basic algorithm against other styles of allocators.
class RABasic : public MachineFunctionPass,
                public RegAllocBase,
                private LiveRangeEdit::Delegate {
  // context
  MachineFunction *MF;

  // state
  std::unique_ptr<Spiller> SpillerInstance;
  std::priority_queue<LiveInterval*, std::vector<LiveInterval*>,
                      CompSpillWeight> Queue;

  // Scratch space.  Allocated here to avoid repeated malloc calls in
  // selectOrSplit().
  BitVector UsableRegs;

  bool LRE_CanEraseVirtReg(unsigned) override;
  void LRE_WillShrinkVirtReg(unsigned) override;
  void printInterGraph(DenseMap <unsigned, std::vector<unsigned>> *InterGraph);
  void checkInterferences(DenseMap <unsigned, std::vector<unsigned>> *InterGraph, LiveInterval *VirtReg);
  bool isSEO(SmallVector<unsigned, 8> SEO, DenseMap <unsigned, std::vector<unsigned>> *InterGraph);
  bool MCS(DenseMap <unsigned, std::vector<unsigned>> *InterGraph);


public:
  RABasic();

  /// Build the Interference Graph of the program
  void buildInterGraph();

  /// Get the number of vertices and edges of a graph
  std::pair<unsigned,unsigned> getInterGraphSize(DenseMap <unsigned, std::vector<unsigned>> *InterGraph);

  /// Return the pass name.
  StringRef getPassName() const override { return "Basic Register Allocator"; }

  /// RABasic analysis usage.
  void getAnalysisUsage(AnalysisUsage &AU) const override;

  void releaseMemory() override;

  Spiller &spiller() override { return *SpillerInstance; }

  void enqueue(LiveInterval *LI) override {
    Queue.push(LI);
  }

  LiveInterval *dequeue() override {
    if (Queue.empty())
      return nullptr;
    LiveInterval *LI = Queue.top();
    Queue.pop();
    return LI;
  }

  unsigned selectOrSplit(LiveInterval &VirtReg,
                         SmallVectorImpl<unsigned> &SplitVRegs) override;

  /// Perform register allocation.
  bool runOnMachineFunction(MachineFunction &mf) override;

  MachineFunctionProperties getRequiredProperties() const override {
    return MachineFunctionProperties().set(
        MachineFunctionProperties::Property::NoPHIs);
  }

  // Helper for spilling all live virtual registers currently unified under preg
  // that interfere with the most recently queried lvr.  Return true if spilling
  // was successful, and append any new spilled/split intervals to splitLVRs.
  bool spillInterferences(LiveInterval &VirtReg, unsigned PhysReg,
                          SmallVectorImpl<unsigned> &SplitVRegs);

  static char ID;
};

char RABasic::ID = 0;

} // end anonymous namespace

char &llvm::RABasicID = RABasic::ID;

INITIALIZE_PASS_BEGIN(RABasic, "regallocbasic", "Basic Register Allocator",
                      false, false)
INITIALIZE_PASS_DEPENDENCY(LiveDebugVariables)
INITIALIZE_PASS_DEPENDENCY(SlotIndexes)
INITIALIZE_PASS_DEPENDENCY(LiveIntervals)
INITIALIZE_PASS_DEPENDENCY(RegisterCoalescer)
INITIALIZE_PASS_DEPENDENCY(MachineScheduler)
INITIALIZE_PASS_DEPENDENCY(LiveStacks)
INITIALIZE_PASS_DEPENDENCY(MachineDominatorTree)
INITIALIZE_PASS_DEPENDENCY(MachineLoopInfo)
INITIALIZE_PASS_DEPENDENCY(VirtRegMap)
INITIALIZE_PASS_DEPENDENCY(LiveRegMatrix)
INITIALIZE_PASS_END(RABasic, "regallocbasic", "Basic Register Allocator", false,
                    false)


//Check if a given vertices order is indeed a simplicial elimination order
bool RABasic::isSEO(SmallVector<unsigned, 8> SEO, DenseMap <unsigned, std::vector<unsigned>> *InterGraph){
  //For each vertex v in the SEO, we must check if v is simplicial in the subgraph induced by {v1,...,v}
  std::set<unsigned> Subgraph;
  int Max_Clique_Size = 0;
  int Clique_Size = 0;
  for(auto& v : SEO){
    //Subgraph induced by all vertices previous v plus v itself
    Subgraph.insert(v);
    for(auto& u : Subgraph){
      //We need to check if the vertices u in the graph, (except v itself) form a clique.
      //For each vertex in the subgraph, go through its neighbourhood.
      for(auto& n1 : (*InterGraph)[u]){
        if(Subgraph.find(n1) != Subgraph.end()){
          for(auto& n2 : (*InterGraph)[u]){
            //We don't have to consider the v node itself
            if(n1==n2 || n1==v || n2==v)
              continue;
            if(Subgraph.find(n2) != Subgraph.end()){
              //Search fot n2 in the neighbourhood of n1.
              if(std::find((*InterGraph)[n1].begin(), (*InterGraph)[n1].end(), n2) == (*InterGraph)[n1].end()){
                //If two neighbours of a node are not adjacent, then, there's no clique.
                //Thus v is not simplicial in the subgraph, which implies the simplicial elimination ordering
                //is not valid. Thereafter the graph is not chordal.
                return false;
              }
              else
                Clique_Size = (*InterGraph)[n1].size();
            }
          }
        }
      }
    }
    if(Clique_Size > Max_Clique_Size)
      Max_Clique_Size = Clique_Size;
    Clique_Size = 0;
  }
  errs() << Max_Clique_Size << " ";
  return true;
}


//Maximum Cardinality Search algorithm to determine of a graph has a simplicial elimination order,
//therefore, that graph is chordal.
bool RABasic::MCS(DenseMap <unsigned, std::vector<unsigned>> *InterGraph){
  DenseMap <unsigned, int> Weights(InterGraph->size());
  DenseMap <unsigned, bool> Visiteds(InterGraph->size());
  SmallVector<unsigned, 8> SEO;

  //Initialize the degree of all vertices with zero and set the vertex as unvisited
  for(auto& vertex: *InterGraph){
    Weights.insert(std::make_pair(vertex.first, 0));
    Visiteds.insert(std::make_pair(vertex.first, false));
  }

  for(auto& vertex: *InterGraph){
    //If we already visited a vertex, there's no need to visit it again
    if(Visiteds[vertex.first])
      continue;
    int v_weight=Weights[vertex.first];
    //Check if this vertex has degree greater or equal than all the other vertices in the graph
    if(std::all_of(Weights.begin(), Weights.end(), [&v_weight] (std::pair<unsigned,int> x) {return x.second<=v_weight;})){
      //Add this vertex in the order
      SEO.push_back(vertex.first);
      //For all of its neighbours, increase degree by one and mark this vertex as visited
      for(auto& neighbour: vertex.second){
        if(Visiteds[neighbour])  //If a vertex belongs to V(G) inter N(v)
          continue;
        Weights[neighbour] += 1;
        Visiteds[vertex.first] = true;
      }
    }
  }

  //In order to state that a graph is chordal, we need to check if the produced order is indeed a
  //simplicial elimination ordering
  return isSEO(SEO, InterGraph);
}

//Print the Interference Graph in a adjacency list representation
void RABasic::printInterGraph(DenseMap <unsigned, std::vector<unsigned>> *InterGraph){
  errs() <<"Interference Graph:\n";
  for(auto&& [vertex,neighbours]: *InterGraph){
    errs() << printReg(vertex, TRI) << "-> ";
    for(unsigned &neighbour: neighbours)
      errs() << printReg(neighbour, TRI) << " ";
     errs() <<"\n";
  }
}

//Get the number of vertices and number of edges of a Interference Graph
std::pair<unsigned,unsigned> RABasic::getInterGraphSize(DenseMap <unsigned, std::vector<unsigned>> *InterGraph){
  unsigned numVertices = InterGraph->size();
  unsigned numEdges = 0;
  for(auto&& [vertex,neighbours]: *InterGraph)
    numEdges+=neighbours.size();

  numEdges/=2;
  return std::make_pair(numVertices, numEdges);
}

//Check for interference between a virtual register and all the others
void RABasic::checkInterferences(DenseMap <unsigned, std::vector<unsigned>> *InterGraph, LiveInterval *VirtReg){
  //For all the other virtual registers
  for(unsigned i = 0, e = MRI->getNumVirtRegs(); i != e; ++i) {
    unsigned AnotherReg = TargetRegisterInfo::index2VirtReg(i);
    if(MRI->reg_nodbg_empty(AnotherReg))
      continue;
    LiveInterval *AnotherVirtReg = &LIS->getInterval(AnotherReg);

    //If we are handling the same register, skip it
    if(AnotherVirtReg->reg == VirtReg->reg)
      continue;

    if(AnotherVirtReg->overlaps(*VirtReg)){
      //If two regs interfere, create an edge in the graph
      (*InterGraph)[VirtReg->reg].push_back(AnotherVirtReg->reg);
    }
  }
}

//Main function for building the Interference Graph
void RABasic::buildInterGraph(){
  DenseMap <unsigned, std::vector<unsigned>> InterGraph;
  //errs() <<"Building interference graph...\n";
  std::vector<unsigned> Edges;
	for(unsigned i = 0, e = MRI->getNumVirtRegs(); i != e; ++i) {
  	//Register ID
  	unsigned Reg = TargetRegisterInfo::index2VirtReg(i);

  	//If is not a register used/defined only in debug instructions
  	if(MRI->reg_nodbg_empty(Reg))
      continue;

    //Create a node for this virtual Register
    InterGraph.insert(std::make_pair(Reg,Edges));

    //Get the respective LiveInterval
    LiveInterval *VirtReg = &LIS->getInterval(Reg);

    //Check for interferences among the virtual registers
    checkInterferences(&InterGraph, VirtReg);

  }
  printInterGraph(&InterGraph);
  MCS(&InterGraph) ? errs() << "\nGraph Chordal\n" : errs() << "\nGraph Non-Chordal\n";
  std::pair<unsigned,unsigned> IGS=getInterGraphSize(&InterGraph);
  errs() << "Number of Vertices: "<< IGS.first <<" \nNumber of Edges: " << IGS.second <<"\n";
}


bool RABasic::LRE_CanEraseVirtReg(unsigned VirtReg) {
  LiveInterval &LI = LIS->getInterval(VirtReg);
  if (VRM->hasPhys(VirtReg)) {
    Matrix->unassign(LI);
    aboutToRemoveInterval(LI);
    return true;
  }
  // Unassigned virtreg is probably in the priority queue.
  // RegAllocBase will erase it after dequeueing.
  // Nonetheless, clear the live-range so that the debug
  // dump will show the right state for that VirtReg.
  LI.clear();
  return false;
}

void RABasic::LRE_WillShrinkVirtReg(unsigned VirtReg) {
  if (!VRM->hasPhys(VirtReg))
    return;

  // Register is assigned, put it back on the queue for reassignment.
  LiveInterval &LI = LIS->getInterval(VirtReg);
  Matrix->unassign(LI);
  enqueue(&LI);
}

RABasic::RABasic(): MachineFunctionPass(ID) {
}

void RABasic::getAnalysisUsage(AnalysisUsage &AU) const {
  AU.setPreservesCFG();
  AU.addRequired<AAResultsWrapperPass>();
  AU.addPreserved<AAResultsWrapperPass>();
  AU.addRequired<LiveIntervals>();
  AU.addPreserved<LiveIntervals>();
  AU.addPreserved<SlotIndexes>();
  AU.addRequired<LiveDebugVariables>();
  AU.addPreserved<LiveDebugVariables>();
  AU.addRequired<LiveStacks>();
  AU.addPreserved<LiveStacks>();
  AU.addRequired<MachineBlockFrequencyInfo>();
  AU.addPreserved<MachineBlockFrequencyInfo>();
  AU.addRequiredID(MachineDominatorsID);
  AU.addPreservedID(MachineDominatorsID);
  AU.addRequired<MachineLoopInfo>();
  AU.addPreserved<MachineLoopInfo>();
  AU.addRequired<VirtRegMap>();
  AU.addPreserved<VirtRegMap>();
  AU.addRequired<LiveRegMatrix>();
  AU.addPreserved<LiveRegMatrix>();
  MachineFunctionPass::getAnalysisUsage(AU);
}

void RABasic::releaseMemory() {
  SpillerInstance.reset();
}


// Spill or split all live virtual registers currently unified under PhysReg
// that interfere with VirtReg. The newly spilled or split live intervals are
// returned by appending them to SplitVRegs.
bool RABasic::spillInterferences(LiveInterval &VirtReg, unsigned PhysReg,
                                 SmallVectorImpl<unsigned> &SplitVRegs) {
  // Record each interference and determine if all are spillable before mutating
  // either the union or live intervals.
  SmallVector<LiveInterval*, 8> Intfs;

  // Collect interferences assigned to any alias of the physical register.
  for (MCRegUnitIterator Units(PhysReg, TRI); Units.isValid(); ++Units) {
    LiveIntervalUnion::Query &Q = Matrix->query(VirtReg, *Units);
    Q.collectInterferingVRegs();
    for (unsigned i = Q.interferingVRegs().size(); i; --i) {
      LiveInterval *Intf = Q.interferingVRegs()[i - 1];
      if (!Intf->isSpillable() || Intf->weight > VirtReg.weight)
        return false;
      Intfs.push_back(Intf);
    }
  }
  LLVM_DEBUG(dbgs() << "spilling " << printReg(PhysReg, TRI)
                    << " interferences with " << VirtReg << "\n");
  assert(!Intfs.empty() && "expected interference");

  // Spill each interfering vreg allocated to PhysReg or an alias.
  for (unsigned i = 0, e = Intfs.size(); i != e; ++i) {
    LiveInterval &Spill = *Intfs[i];

    // Skip duplicates.
    if (!VRM->hasPhys(Spill.reg))
      continue;

    // Deallocate the interfering vreg by removing it from the union.
    // A LiveInterval instance may not be in a union during modification!
    Matrix->unassign(Spill);

    // Spill the extracted interval.
    LiveRangeEdit LRE(&Spill, SplitVRegs, *MF, *LIS, VRM, this, &DeadRemats);
    spiller().spill(LRE);
  }
  return true;
}

// Driver for the register assignment and splitting heuristics.
// Manages iteration over the LiveIntervalUnions.
//
// This is a minimal implementation of register assignment and splitting that
// spills whenever we run out of registers.
//
// selectOrSplit can only be called once per live virtual register. We then do a
// single interference test for each register the correct class until we find an
// available register. So, the number of interference tests in the worst case is
// |vregs| * |machineregs|. And since the number of interference tests is
// minimal, there is no value in caching them outside the scope of
// selectOrSplit().
unsigned RABasic::selectOrSplit(LiveInterval &VirtReg,
                                SmallVectorImpl<unsigned> &SplitVRegs) {
  // Populate a list of physical register spill candidates.
  SmallVector<unsigned, 8> PhysRegSpillCands;

  // Check for an available register in this class.
  AllocationOrder Order(VirtReg.reg, *VRM, RegClassInfo, Matrix);
  while (unsigned PhysReg = Order.next()) {
    // Check for interference in PhysReg
    switch (Matrix->checkInterference(VirtReg, PhysReg)) {
    case LiveRegMatrix::IK_Free:
      // PhysReg is available, allocate it.
      return PhysReg;

    case LiveRegMatrix::IK_VirtReg:
      // Only virtual registers in the way, we may be able to spill them.
      PhysRegSpillCands.push_back(PhysReg);
      continue;

    default:
      // RegMask or RegUnit interference.
      continue;
    }
  }

  // Try to spill another interfering reg with less spill weight.
  for (SmallVectorImpl<unsigned>::iterator PhysRegI = PhysRegSpillCands.begin(),
       PhysRegE = PhysRegSpillCands.end(); PhysRegI != PhysRegE; ++PhysRegI) {
    if (!spillInterferences(VirtReg, *PhysRegI, SplitVRegs))
      continue;

    assert(!Matrix->checkInterference(VirtReg, *PhysRegI) &&
           "Interference after spill.");
    // Tell the caller to allocate to this newly freed physical register.
    return *PhysRegI;
  }

  // No other spill candidates were found, so spill the current VirtReg.
  LLVM_DEBUG(dbgs() << "spilling: " << VirtReg << '\n');
  if (!VirtReg.isSpillable())
    return ~0u;
  LiveRangeEdit LRE(&VirtReg, SplitVRegs, *MF, *LIS, VRM, this, &DeadRemats);
  spiller().spill(LRE);

  // The live virtual register requesting allocation was spilled, so tell
  // the caller not to allocate anything during this round.
  return 0;
}

bool RABasic::runOnMachineFunction(MachineFunction &mf) {
  LLVM_DEBUG(dbgs() << "********** BASIC REGISTER ALLOCATION **********\n"
                    << "********** Function: " << mf.getName() << '\n');

  MF = &mf;

  RegAllocBase::init(getAnalysis<VirtRegMap>(),
                     getAnalysis<LiveIntervals>(),
                     getAnalysis<LiveRegMatrix>());

  calculateSpillWeightsAndHints(*LIS, *MF, VRM,
                                getAnalysis<MachineLoopInfo>(),
                                getAnalysis<MachineBlockFrequencyInfo>());

  buildInterGraph();


  SpillerInstance.reset(createInlineSpiller(*this, *MF, *VRM));

  allocatePhysRegs();
  postOptimization();

  // Diagnostic output before rewriting
  LLVM_DEBUG(dbgs() << "Post alloc VirtRegMap:\n" << *VRM << "\n");

  releaseMemory();
  return true;
}

FunctionPass* llvm::createBasicRegisterAllocator()
{
  return new RABasic();
}
