# Read mapping with insertions
CS CM124 Project #3   
Byron Lutz   
May 27, 2014   

## Naive implementation
Works most of the time but has a few classes of bugs:

* An insert is considered ended once it encounters a matching base. Even if the insert has not actually ended
* Only works for inserts up to length 5
