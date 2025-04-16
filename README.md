STOP: Simulated TCSPC On Proteins
=================================

do you like to shoot lasers at proteins? have you ever thought that maybe the whole process of getting hold of a protein and a laser in real life to do that was a bit too much effort? well, do i have the software for you! it simulates the process of firing an expensive laser at your weird little protein, and you can do it from the comfort of your own desk or couch or wherever you like, you wastrel.

on a more serious note, this is some fortran code to simulate TCSPC experiments on arbitrary proteins (with some constraints, which i'll go into). the idea is that you set up your protein using JSON or python, there's some python code to write all the relevant parameters to fortran-friendly files and then an MPI fortran kernel to do the heavy lifting. once the simulated experiments are done, reconvolution fits are performed to get amplitude-weighted lifetimes and so on.

FAQS
====

Q: where is the test suite, the continuous integration, all that kind of stuff?  
A: what are you, a coward?
