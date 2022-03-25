# mIcal track reconstruction

The code `ino_digi_read1.C` is to reconstruct tracks in mIcal using `TMinuit`.

The function `PropagateTrack` propagates a muon with the given input parameters.

There are two method to prapagate. Defined with the flag `PropagateTrack`.

If the flag `isSimData` is defined, then it works for the simulated events.

Structure of events for data and simulation are given in the `RPCEve.h` and `T2.h` class files, respectively.

Muon energy loss file: `muon_energy_loss_in_fe.txt`

The detector alignment data is in the ASCII files, i.e. `mIcal_correction_PosTime_20210702.txt  mIcal_correction_sim.txt  mIcal_correction_TDCstrp_20210702.txt`

Input file location is hard coded in the code.
```
For Simulated: `/var/nfscondor/surya/sim/`
For Data: `/var/nfscondor/surya/`
```

Output file location is hard coded in the code: `./temp/`.

One might change it in the code, if the input/output file location differs.

A bunch of input files may be present in the above directory. One might use them directly.  

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
./ino_digi_read1 BRPCv4t_evtraw_20181217_162121.root 0 100 0
./ino_digi_read1 corsika76300_FLUKA_SIBYLL_3dFlux_20220105av_trg5of8_20220105aa_25Cr_digi.root 0 100 0
```
