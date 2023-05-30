import os

import pymongo
import urllib.parse
from pymongo.server_api import ServerApi
from dotenv import load_dotenv
from tqdm import tqdm

load_dotenv()

USERNAME = urllib.parse.quote_plus(os.environ.get("MONGO_USERNAME"))
PASSWORD = urllib.parse.quote_plus(os.environ.get("MONGO_PASSWORD"))
connection_url = f"mongodb+srv://{USERNAME}:{PASSWORD}@socialmediadatabase.gihvf.mongodb.net/?retryWrites=true&w=majority"


client = pymongo.MongoClient(f"mongodb+srv://{USERNAME}:{PASSWORD}@socialmediadatabase.gihvf.mongodb.net/SocialMediaDatabase?retryWrites=true&w=majority", server_api=ServerApi('1'))
# client = pymongo.MongoClient(f"mongodb+srv://socialmediadatabase.gihvf.mongodb.net/SocialMediaDatabase?retryWrites=true&w=majority", server_api=ServerApi('1'))


DATABASE_NAME = "SocialMediaCaregivingResearch"
COLLECTION_NAME = "AlzConnected"
db = client[DATABASE_NAME]
collection = db[COLLECTION_NAME]


def parse_document(doc:dict, collection:str):
    return {
        "id": str(doc["_id"]),
        "date":doc["date"],
        "collection":collection,
        "url":doc["url"],
        "title":doc["title"],
        "body":doc["body"]
    }


data = {}

for doc in tqdm(collection.find({}),desc="load collection from db"):
    doc = parse_document(doc,COLLECTION_NAME)
    url = doc["url"]
    if(url not in data):
        data[url] = []
    data[url].append(doc)

sorted_data = {url:sorted(docs,key=lambda x:x["date"]) for url,docs in data.items()}


import pickle

with open(f"./{COLLECTION_NAME}.pkl","wb") as output_file:
    pickle.dump(sorted_data,output_file)


