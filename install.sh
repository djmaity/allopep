###################
# Install APOP    #
###################
git clone https://github.com/Ambuj-UF/APOP.git

###################
# Install GAPS    #
###################
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

###################
# Install PepGLAD #
###################
# Clone the PepGLAD GitHub repo
git clone https://github.com/THUNLP-MT/PepGLAD.git

# Download the PepGLAD weights
wget https://github.com/THUNLP-MT/PepGLAD/releases/download/v1.0/checkpoints.zip

# Extract the PepGLAD weights into PepGLAD directory
unzip checkpoints.zip "checkpoints/*" -d PepGLAD
