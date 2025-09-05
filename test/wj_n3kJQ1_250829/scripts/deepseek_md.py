import os
import requests
from datetime import datetime
import pandas as pd

DEEPSEEK_API_URL = "https://api.deepseek.com/v1/chat/completions"
DEEPSEEK_API_KEY = "sk-e059549291fd4ba891a0690d110257d0"  # Replace with your actual API key

def generate_markdown_with_deepseek():
    # Find files in current directory
    txt_files = sorted(f for f in os.listdir() if f.endswith('.txt'))
    tsv_files = sorted(f for f in os.listdir() if f.endswith('.tsv'))
    png_files = sorted(f for f in os.listdir() if f.endswith('.png'))
    
    # Prepare content for DeepSeek API
    content_parts = []
    
    # Process text files
    if txt_files:
        content_parts.append("## Text Content\n\n")
        for txt_file in txt_files:
            with open(txt_file, 'r') as f:
                content_parts.append(f"### {txt_file}\n```\n{f.read()}\n```\n\n")
    
    # Process TSV files as tables
    if tsv_files:
        content_parts.append("## Data Tables\n\n")
        for tsv_file in tsv_files:
            df = pd.read_csv(tsv_file, sep='\t')
            content_parts.append(f"### {tsv_file}\n\n{df.to_markdown(index=False)}\n\n")
    
    # Process PNG files
    if png_files:
        content_parts.append("## Visualizations\n\n")
        for i, png_file in enumerate(png_files, 1):
            content_parts.append(f"![Figure {i}]({png_file})\n*Figure {i}: {os.path.splitext(png_file)[0]}*\n\n")
    
    # Combine all content
    full_content = "# Project Documentation\n\n" + "".join(content_parts)
    full_content += f"Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
    
    # Call DeepSeek API to generate final markdown
    headers = {
        "Authorization": f"Bearer {DEEPSEEK_API_KEY}",
        "Content-Type": "application/json"
    }
    
    payload = {
        "model": "deepseek-chat",
        "messages": [
            {
                "role": "system",
                "content": "You are a helpful assistant that generates well-formatted markdown documentation."
            },
            {
                "role": "user",
                "content": f"Please format this content as a professional markdown document:\n\n{full_content}"
            }
        ],
        "temperature": 0.7
    }
    
    response = requests.post(DEEPSEEK_API_URL, json=payload, headers=headers)
    response.raise_for_status()
    
    # Save the generated markdown
    with open("project_documentation.md", "w") as f:
        f.write(response.json()["choices"][0]["message"]["content"])
    
    print("Markdown documentation generated: project_documentation.md")

if __name__ == "__main__":
    generate_markdown_with_deepseek()