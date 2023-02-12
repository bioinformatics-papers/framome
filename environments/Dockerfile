FROM continuumio/miniconda3:22.11.1

# Create environtment
COPY environment.yml /setup/environment.yml
RUN conda env create -f /setup/environment.yml
RUN echo "source activate framome_paper" > ~/.bashrc
ENV PATH /opt/conda/envs/myenv/bin:$PATH

# Execute
CMD [ "bash" ]

