# Use Ubuntu as the base image
FROM ubuntu:latest

# Install C++, CMake, Eigen, and OpenMP
RUN apt-get update && apt-get install -y \
    g++ \
    cmake \
    libeigen3-dev \
    libomp-dev

# Set the working directory inside the container
WORKDIR /usr/src/app

# Copy the CMakeLists.txt, source code, and the include directory into the container
COPY CMakeLists.txt .
COPY src/ ./src
COPY include/ ./include
COPY bindings/ ./bindings
COPY external/ ./external

# Build the project
# Assuming your main executable is defined in your CMakeLists.txt
RUN cmake . && make

# Run the compiled binary
# Replace `your_executable_name` with the actual name of your binary
CMD ["./clipperplus_test"]
