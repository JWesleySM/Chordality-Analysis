# Chordality Analysis

This repository contains source code to build a interference graph of a program. This implementation will also give you the size of the graph, amount of vertices and edges, the size of a maximum clique and information this graph is chordal or not. All the logic was inserted into the LLVM Basic Register Allocator. The following methods are added:

```
buildInterGraph() -> build the interference graph of a program.
printInterGraph() -> print the interference graph.
checkInterferences() -> check if a virtual register is simultaneously alive as any other.
MCS() -> Maximum Cardinality Search algorithm to create a Simplicial Elimination Ordering for a graph
isSEO() -> check if a sequence is indeed a simplicial elimination ordering. It also constructs the 
cliques and prints the size of the largest one.
getInterGraphSize() -> get the number of vertices and edges in a graph.
```

In order to use it, [download and build](http://releases.llvm.org) on your system. Then replace the file located at <your-path-to-llvm/llvm/lib/CodeGen/RegAllocBasic.cpp for this one and recompile the [LLVM static compiler](https://llvm.org/docs/CommandGuide/llc.html).

```
cd <path-to-your-llvm>
cd build
make llc
```
To use it on a program, first compile the source code (.c file) with Clang, and call the static compiler llc on the bytecode. For example

```
add here information to call llc on a C program
```

This implementation is part of the [Angha Project](http://cuda.dcc.ufmg.br/angha/home). By following this link, you can find a very nice analysis on chordal graphs in the tab "Analyses". 
