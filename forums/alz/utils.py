def parse_document(document: dict, collection: str):
    doc = {
        "id": str(document["_id"]),
        "collection": collection,
        "url": document["url"],
        "title": document["title"],
        "body": document["body"],
    }
    if collection == "AlzConnected":
        doc["date"] = document["date"]
    elif collection == "AlzConnected-v2":
        doc["date"] = document["datetime"]
    return doc
