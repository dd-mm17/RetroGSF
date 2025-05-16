![Project Logo](assets/banner.png)

![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
RetroGSF
</h1>

<br>


Predict and evaluates a solvent for a retrosynthesized reaction

## 🔥 Usage

```python
from mypackage import main_func

# One line to rule them all
result = main_func(data)
```
To use the applet...


## 👩‍💻 Installation

For the first time you will need to create a new environment, you may also give the environment a different name. 

```
conda create -n retrogsf python=3.10 
```

Then, activate your environment to install the package using pip:
```
conda activate retrogsf
(conda_env) $ pip install .
```

## Additional installations
In order to use the streamlit applet you will need to have an Aizynthfinder config.yml file along with a Google AI Studios API keys.

### config.yml file
The config.yml file can be created using (where ```my_folder``` is the folder that you want download to): 
```
download_public_data my_folder
```

More information can be found found here: [Documentation](https://molecularai.github.io/aizynthfinder/#) or here [GithHub](https://github.com/MolecularAI/aizynthfinder?tab=readme-ov-file)

The path to your config.yml file will need to be updated in the code at the following areas:


### Google API key

1.) Open the following link : (https://aistudio.google.com/app/apikey)

2.) Click on "Get API key"

3.) Create a .env file in the root folder with the following text : ```GEMINI_API_KEY="YOUR_API_KEY"```


## 🛠️ Development installation

Initialize Git (only for the first time). 

Note: You should have create an empty repository on `https://github.com:dd-mm17/retrogsf`.

```
git init
git add * 
git add .*
git commit -m "Initial commit" 
git branch -M main
git remote add origin git@github.com:dd-mm17/retrogsf.git 
git push -u origin main
```

Then add and commit changes as usual. 

To install the package, run

```
(retrogsf) $ pip install -e ".[test,doc]"
```

### Run tests and coverage

```
(conda_env) $ pip install tox
(conda_env) $ tox
```



