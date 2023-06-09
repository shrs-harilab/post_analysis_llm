import os
import pymongo
import urllib.parse
from pymongo.server_api import ServerApi
from dotenv import load_dotenv
from tqdm import tqdm
import pickle
from utils import parse_document


load_dotenv()

USERNAME = urllib.parse.quote_plus(os.environ.get("MONGO_USERNAME"))
PASSWORD = urllib.parse.quote_plus(os.environ.get("MONGO_PASSWORD"))

client = pymongo.MongoClient(
    f"mongodb+srv://{USERNAME}:{PASSWORD}@socialmediadatabase.gihvf.mongodb.net/SocialMediaDatabase?retryWrites=true&w=majority",
    server_api=ServerApi("1"),
)

DATABASE_NAME = "SocialMediaCaregivingResearch"
main = "AlzConnected"
COLLECTION_NAMES = [main, "AlzConnected-v2"]
db = client[DATABASE_NAME]
collections = [db[name] for name in COLLECTION_NAMES]

data = {}

for collection in collections:
    for doc in tqdm(collection.find({}), desc=f"loading {collection.name} from db"):
        doc = parse_document(doc, collection.name)
        url = doc["url"]
        if url not in data:
            data[url] = []
        data[url].append(doc)

sorted_data = {url: sorted(docs, key=lambda x: x["date"]) for url, docs in data.items()}


with open(os.path.join(os.path.dirname(__file__), f"{main}.pkl"), "wb") as output_file:
    pickle.dump(sorted_data, output_file)
