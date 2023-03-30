#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --job-name=flow-test-run

echo "This should run your application."
./flow test_images/500x500_0.bmp test_images/500x500_1.bmp
