#!/bin/bash

set -e

while getopts "sb:" opt; do
    case $opt in
        s) use_sra=1 ;;
        b) branch="$OPTARG" ;;
        *) echo "Usage: $0 [-s] [-b <branch_name>]" && exit 1
    esac
done
shift $(($OPTIND - 1))

if [ "$branch" == "" ] ; then
    branch="master"
fi

set -x
yum install -y git zip unzip pandoc

git clone https://github.com/BenLangmead/bowtie.git
if [ $? -ne 0 ] ; then
    echo "Unable to clone bowtie repo"
    exit 1
fi

cd bowtie

pwd && ls
git branch -a | grep "$branch" 2>&1 > /dev/null
if [ $? -ne 0 ] ; then
    echo "branch '$branch' does not exist"
    exit 1
else
     git checkout "$branch"
fi

if [ $use_sra -eq 1 ] ; then
    # this variant is needed to compile ncbi-vdb
    source /hbb/activate
    make sra-deps
fi

# this variant creates static binaries with PIC
source /hbb_exe_gc_hardened/activate

make bowtie-bin.zip
if [ $? -ne 0 ] ; then
    echo "Unable to create bowtie package"
    exit 1
fi
echo "Running libcheck..."
libcheck bowtie-{align,build,inspect}-*

echo "Running hardening check..."
hardening-check -b bowtie-{align,build,inspect}-*

echo "Copying binary package"
cp /bowtie/*.zip /io
