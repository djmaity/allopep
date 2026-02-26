git clone https://github.com/Ambuj-UF/APOP.git

git clone https://github.com/google-deepmind/alphafold3.git

git clone https://github.com/hongliangduan/GAPS.git
cd GAPS
conda create -n gaps python=3.10
conda activate gaps
conda install pytorch==1.13.1 torchvision==0.14.1 torchaudio==0.13.1 pytorch-cuda=11.7 -c pytorch -c nvidia
conda install pandas scikit-learn tqdm h5py gemmi

