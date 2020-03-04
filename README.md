# Chordality Analysis

This repository contains source code to build a interference graph of a program. This implementation will also give you the size of the graph, amount of vertices and edges and information this graph is chordal or not. All the logic was inserted into the LLVM Basic Register Allocator. The following methods are added:

```
buildInterGraph() -> build the interference graph of a program.
printInterGraph() -> print the interference graph.
checkInterferences() -> check if a virtual register is simultaneously alive as any other.
MCS() -> Maximum Cardinality Search algorithm to create a Simplicial Elimination Ordering for a graph
isSEO() -> check if a sequence is indeed a simplicial elimination ordering.
getInterGraphSize() -> get the number of vertices and edges in a graph.
```

In order to use it, [download and build](http://releases.llvm.org) on your system. Then replace the file located at <your-path-to-llvm/llvm/lib/CodeGen/RegAllocBasic.cpp for this one and recompile the [LLVM static compiler](https://llvm.org/docs/CommandGuide/llc.html).

```
cd <path-to-your-llvm>
cd build
make llc
```
To use it on a program, first compile the source code (.c file) with Clang, and call the static compiler llc on the bytecode. For example:

```
clang -emit-llvm program.c -c -o program.bc
llc -regalloc=basic program.bc -o program.s
```
It is worth remembering that this code handle a LLVM Machine Function. Then the output contains the name of the virtual register, not the variables in the source code. After running llc, you will see an output like this:
```
Interference Graph:
%26-> 
%7-> %6 %8 %9 %10 %11 %12 %13 
%21-> 
%9-> %6 %7 %8 %10 %11 %12 %13 
%11-> %6 %7 %8 %9 %10 %12 %13 
%6-> %7 %8 %9 %10 %11 %12 %13 
%13-> %6 %7 %8 %9 %10 %11 
%8-> %6 %7 %9 %10 %11 %12 %13 
%34-> 
%15-> 
%29-> 
%10-> %6 %7 %8 %9 %11 %12 %13 
%17-> 
%24-> 
%12-> %6 %7 %8 %9 %10 %11 
%19-> 

Graph Chordal
Number of Vertices: 16 
Number of Edges: 27
```
*Obs: the graph here is represented as an adjacency list. Thus it is read this way:*

*vertex-> its neighbours*

If for any reason you don't want some of the output information, comment the corresponding lines 307 to 309 at the **RegAllocBasic.cpp** file.

This implementation is part of the [Angha Project](http://cuda.dcc.ufmg.br/angha/home). By following this link, you can find a very nice analysis on chordal graphs in the tab "Analyses". 
