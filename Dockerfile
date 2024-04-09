FROM --platform=linux/amd64 python:3.10.9-slim AS build

MAINTAINER Michael Sorenson <masorens@ucsd.edu>

WORKDIR /usr/src/app

COPY requirements.txt ./

RUN pip install --upgrade pip
RUN pip install -r requirements.txt

COPY . .

EXPOSE 80

CMD ["python", "app.py", "80"]