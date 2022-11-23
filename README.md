# mIcal track reconstruction

The code `ino_digi_read1.C` is to reconstruct tracks in mIcal using `TMinuit`.

The function `PropagateTrack` propagates a muon with the given input parameters.

There are two method to prapagate. Defined with the flag `PropagateTrack`.

If the flag `isSimData` is defined, then it works for the simulated events.

Structure of events for data and simulation are given in the `RPCEve.h` and `T2.h` class files, respectively.

Muon energy loss file: `muon_energy_loss_in_fe.txt`

The detector alignment data is in the ASCII files, i.e. `mIcal_correction_PosTime_20210702.txt  mIcal_correction_sim.txt  mIcal_correction_TDCstrp_20210702.txt`

Input file location is hard coded in the code `./input_files/`.

Output file location is hard coded in the code: `./temp/`.

One might change it in the code, if the input/output file location differs.

A bunch of input files may be present in `/var/nfscondor/surya/` and `/var/nfscondor/surya/sim/`. One might use them.  

Arguments to the executable:
```
0. executable itself
1. input file
2. starting event number (starts from 0)
3. end event number
4. output file file suffix (any positive integer upto 99999); useful for batch submission.
```

Source in sim01: `source env.sh`

Compile: `./execute`

Run:
```
For simulated events: ./ino_digi_read1 test_digi.root 0 100 0
For Data:             ./ino_digi_read1 BRPCv4t_evtraw_20181226_192148.root 0 100 0
```

Batch Submission: The following scrips are for submitting jobs to htcondor.
```
surya_job_data.jdl surya_job_sim.jdl
```
Change the `initialdir` before submitting jobs.

## Collated file
Check the folder `collated` for details.