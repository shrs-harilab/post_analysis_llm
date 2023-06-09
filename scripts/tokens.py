import tiktoken
import os

if __name__ == "__main__":
    model = "gpt-3.5-turbo"
    encoding = tiktoken.encoding_for_model(model)
    file = os.path.join(os.path.dirname(__file__), "prompt.txt")
    with open(file, "r") as f:
        prompt = f.read()
    encoded = encoding.encode(prompt)

    print(f"Number of words in prompt: {len(prompt.split())}")
    print(f"Number of characters in prompt: {len(prompt)}")
    print(f"Number of tokens: {len(encoded)}")
