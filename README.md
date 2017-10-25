this module allows encoding kmers of length up to 31 into uint64's

it is as fast as possible and allows reverse-complenting on encoded data
and adding a single base to the right end or a complemented base for the
left end as is needed when sliding along a sequence to generate all kmers.
