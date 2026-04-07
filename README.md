Finds an ideal QM region partition within an Atoms object, for QM/MM models.

Ideal here is defined as having the lowest edge cut, that is the fewest
bonds between the QM region and the surrounding MM region.

An example script is given after the function definitions.

How the algorithm works
-----------------------
First, a trial QM region is selected by naively picking neighbours around a
seed selection until n atoms are selected. Then, atoms selected as QM are
swapped with MM atoms repeatedly. At each step, the swaps which result in the
lowest edge cuts are determined, and then one is taken randomly. This process
is repeatedly for a number of steps and then the best results are collected.
This is repeated across several threads. Finally, the results from each thread
are compared, and the best solution saved for analysis.

Further Work
------------
- Have criteria to prioritise QM regions over others with equal edge cuts.
ie. average step count from seed.
- Allow users to add atoms that are always excluded from the QM region.
- Hypothetically - bond weighting.
