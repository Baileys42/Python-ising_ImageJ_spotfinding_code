import os
from pathlib import Path

def rename_files(folder: str, name: str):
    folder = Path(folder)
    files = sorted(folder.iterdir())
    for i, f in enumerate(files, start=1):
        f.rename(folder / f"{name}_{i}{f.suffix}")


folder = 'B:/Chaperone_subgroup/Bailey/2026/260424_488flash_532image_542A8/Fluc488_10nMA8542_3uMA8_2uMA2/Export/Chap'
name = 'chap'

rename_files(folder, name)