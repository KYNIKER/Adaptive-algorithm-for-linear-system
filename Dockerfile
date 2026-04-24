# official Julia runtime as a parent image
FROM julia:1.12.6

# install external dependencies:
# - make & C compiler  (for building CRlibm)
# - QT  (for GR.jl)
# - LaTeX & divpng  (for LaTeXStrings.jl)
RUN apt-get update && apt-get -qy install make gcc libqt5widgets5 texlive-latex-base dvipng

# set working directory
WORKDIR /react

# copy current directory into container
COPY . /react
RUN chmod -R 777 /react
