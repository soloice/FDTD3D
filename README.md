FDTD(Finite-difference time-domain) method is a numerical analysis technique used for modeling computational electrodynamics.  The idea is rather simple, but this method involves a lot of computation, which makes it sometimes intolerably slow to run on a typical PC.  As the development of parallel computation, we may run FDTD program on a distributed system.  Though some business software has good implementation on this problem, it is really difficult to find an open-source implementation.  What's more, the numerical simulation of the spread of electromagnetic field is a pretty basic problem, so that many researchers may want to use it.  Therefore, we decided to work on this problem and publish our code.

Our code turns to have a good scalability.  Experiments show that we achieves an 6X speedup on a cluster consisting of 8 computing nodes.

Author: 
Chong Ruan: OpenMP & MPI parallelization, inter-process communication.
Weiqiang Zhu: serial program implementaion, generating test data and result visualization.

Feel free to report any bugs to ruanchong_ruby@163.com.  We wish you a good luck!
