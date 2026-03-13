git clone https://github.com/Ambuj-UF/APOP.git

#################
# Install GAPS
#################
# Clone the GAPS GitHub repo
git clone https://github.com/hongliangduan/GAPS.git

# Go into the cloned GAPS directory
cd GAPS

# Convert GAPS into a python package
touch __init__.py

# Apply custom code edit to support newer dependency versions
git apply ../patches/GAPS.patch

# Go back to the AlloPep directory
cd ..


# Folowing packages are already being installed by conda environment.yml
# conda create -n gaps python=3.10
# conda activate gaps
# conda install pytorch==1.13.1 torchvision==0.14.1 torchaudio==0.13.1 pytorch-cuda=11.7 -c pytorch -c nvidia
# conda install pandas scikit-learn tqdm h5py gemmi

#################
# Install PepGLAD
#################
# Clone the PepGLAD GitHub repo
git clone https://github.com/THUNLP-MT/PepGLAD.git

# Download the PepGLAD weights
wget https://github.com/THUNLP-MT/PepGLAD/releases/download/v1.0/checkpoints.zip

# Extract the PepGLAD weights into PepGLAD directory
unzip checkpoints.zip "checkpoints/*" -d PepGLAD

# Apply PepGLAD patch to allow PepGLAD to run with newer dependecies
git apply ../patches/GAPS.patch

git clone https://github.com/google-deepmind/alphafold3.git

