#!/usr/bin/env bash

# Clone repo
cd txt2onto2.0

# Activate environment
conda activate txt2onto2

# Import nltk
pip install nltk
# python import nltk 
# python nltk.download('stopwords')
pip install nltk torch scipy
pip install transformers scikit-learn
pip install plotly

# Fetch raw data
wget -q https://raw.githubusercontent.com/vbilsmith/Metadata/refs/heads/main/data/GREIN_data.csv   -O demo/GREIN_data.csv
wget -q https://raw.githubusercontent.com/vbilsmith/Metadata/refs/heads/main/data/metadata_human.json -O demo/metadata_human.json
wget -q https://raw.githubusercontent.com/vbilsmith/Metadata/refs/heads/main/data/metadata_mouse.json -O demo/metadata_mouse.json

# commit
# echo "demo/metadata_mouse.json" >> .gitignore
# echo "demo/metadata_human.json" >> .gitignore
# echo "results/" >> .gitignore
# git add .gitignore
# git commit -m ""
# git push origin main
