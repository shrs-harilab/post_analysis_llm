import typing
import PyPDF2

def merge_line_to_paragraph(lines:typing.List[str]):
    paragraphs = []
    current_paragraph = []
    for line in lines:
        # here we use newline and sentence end punctuation to detect paragraph boundary.
        # It's not 100% accurate but we can test its performance in prompt later
        if line == "":
            if current_paragraph:
                paragraph = " ".join(current_paragraph)
                paragraphs.append(paragraph)
                current_paragraph = []
        elif line.endswith(".") or line.endswith("?") or line.endswith("!"):
            current_paragraph.append(line)
            paragraph = " ".join(current_paragraph)
            paragraphs.append(paragraph)
            current_paragraph = []
        else:
            current_paragraph.append(line)

    if current_paragraph:
        paragraph = " ".join(current_paragraph)
        paragraphs.append(paragraph)
    return paragraphs

def parse_paragraphs_from_pdf(input_filepath:str):
    reader = PyPDF2.PdfReader(input_filepath)
    lines = []
    for page in reader.pages:
        text = page.extract_text()
        for x in text.splitlines():
            lines.append(x.strip())
    return merge_line_to_paragraph(lines)


if __name__ == '__main__':
    # input_filepath = "/Users/zhendongwang/Downloads/UB-CAM_Training-Manual.pdf"
    # input_filepath = "/Users/zhendongwang/Downloads/FAMCAM_TrainingManual.pdf"
    input_filepath = "/Users/zhendongwang/Downloads/Sour Seven Validation Study.pdf"
    input_filepath = "/Users/zhendongwang/Downloads/nejmcp1605501.pdf"
    for paragraph in parse_paragraphs_from_pdf(input_filepath):
        print(paragraph)
        print("-----")








