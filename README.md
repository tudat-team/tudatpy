# tudat-bundle

## setup

1. Install the contained `environment.yaml` file to satisfy dependencies

````
conda env create -f environment.yaml
````

2. Activate the environment installed in step 1

````
conda activate tudat-bundle
````

3. Determine your `CONDA_PREFIX` path

````
echo $CONDA_PREFIX
````

4. Set the following CMake build configuration

````
-DCMAKE_PREFIX_PATH=<CONDA_PREFIX>
-DCMAKE_CXX_STANDARD=14
-DBoost_NO_BOOST_CMAKE=ON
````

----------------------- ------------------------------------
![Tip](.images/tip.png)\ Using CLion? In CLion, the convention to set CMake arguments
        is to add them to `File>Settings>Build, Execution, Deployment>CMake Options`
----------------------------------------------------------------
