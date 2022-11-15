import os
from pathlib import Path

src = Path(__file__).parent / "src"

for file in [file for file in src.glob('**/*') if file.is_file()]:
    os.system(f'g++ src/{file.name} -o bin/{file.name}')