FROM kyanagis/exprorer_msmd:2025.01.25

RUN mamba install "python>=3.10" -y
RUN pip install rdkit libcoffee
RUN pip install mypy pytest isort ipykernel
RUN apt update 
RUN apt install -y openssh-client

