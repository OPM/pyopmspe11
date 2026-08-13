WHR="test_outputs/paper_convergence"
if [ ! -d "test_outputs" ]; then
    mkdir "test_outputs"
fi
if [ -d $WHR ]; then
    rm -rf $WHR
fi
mkdir $WHR
cp -r convergence/. $WHR
cd $WHR
python3 convergence.py
cd ../..
