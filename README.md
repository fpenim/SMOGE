#Stochastic Modelling of Gene Expression

Proteins carry out the great majority of the catalytic and structural work within
an organism. The RNA templates used in their synthesis determines their identity, and
this is dictated by which genes are transcribed. Therefore, gene expression is the
fundamental determinant of an organism’s nature.

The main objective of this thesis was to develop a stochastic computational
model able to simulate the gene expression phenomena to best evaluate its critical
points including the process of splicing. With this aim, an extended research was
performed looking for already described similar models, identifying their strengths and
weaknesses, approaches, languages and algorithms used.

Here we present a model developed in Java, which is an object-oriented
language or object-oriented programming (OOP) implementing the Gillespie’s
algorithm. The model receives as input an array of 19 biological parameters. Although
at first the model was time consuming and intense, it was optimized reducing the
simulation time in approximately 250%.

With this model we wish to take a closer look at the regulation of gene
expression, evaluating it with greater accuracy. For that purpose we varied the values
of the transcription initiation, RNA degradation and protein degradation rate constants.
The results evidenced that the transcription initiation and RNA degradation present the
same level of control towards the influence of pre-mRNA, mRNA and protein numbers;
in terms of protein numbers, their influence is lower when compared with the protein
degradation constant (which has no influence on RNAs numbers). 