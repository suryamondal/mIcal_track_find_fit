#!/bin/bash 


filelist=(
    # "../outputrootfiles/corsika76300_FLUKA_SIBYLL_512Lakh_7dd7f63.root"
    "../outputrootfiles/corsika76300_FLUKA_SIBYLL_606Lakh_c3b9d93.root"
    # "../outputrootfiles/corsika76300_FLUKA_SIBYLL_614Lakh_4b097e5.root"
    # "../outputrootfiles/corsika76300_FLUKA_SIBYLL_146Lakh_2aaa127.root"
    # "../outputrootfiles/corsika76300_FLUKA_SIBYLL_1665Lakh_529aa33.root"
    # "../outputrootfiles/corsika76300_FLUKA_SIBYLL_135Lakh_0c526af.root"
    "../outputrootfiles/corsika76300_FLUKA_SIBYLL_431Lakh_95de138.root"
    "../outputrootfiles/corsika76300_FLUKA_SIBYLL_739Lakh_3929420.root"
    "../outputrootfiles/corsika76300_FLUKA_SIBYLL_1000Lakh_96cc012.root"
    "../../anal_magnet_20191014_sim01/recodata/test_reco_b249b8e.root"
)

fileTags=(
    "d0\ minimize\ sqrt\ wt"
    "xy\ minimize\ sqrt\ wt"
    "xy\ minimize\ no\ wt"
    "xy\ minimize\ exp\ wt"
    "kalman"
    )

wDir="./"

hDir=`pwd`'/'

cd $wDir

fileArgs='./SigmaComparison.py '

len=${#fileTags[@]}
for (( i=0; i<$len; i++ )); do
    fileArgs+=' '${fileTags[$i]}
    fileArgs+=' '$hDir${filelist[$i]}
done

echo $fileArgs
echo $fileArgs | sh

cd $hDir

echo "Script done"
