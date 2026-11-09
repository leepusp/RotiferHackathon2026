# Functions needing attention

- rdae.extract_envelope()
    - Currently overwrites original seqobj
    - Doesn't notate the new borders if we use expand
- rdbh.hmmscan
    - Doesn't ungap sequences
- rdag.add_coordinates_to_aln()
    - 0 based
- rdag.extend_aln()
    - 1 based, uses previous function
    - Swaps length of the final alignment with lengths of the full sequences instead
 

