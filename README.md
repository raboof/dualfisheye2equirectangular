# Dual Fisheye to Equirectangular Projection Mapping

Many '360' camera's, such as the [dokicam](http://dokicam.com/), consist of 2
fish-eye camera's. 

## Why DIY?

Those camera's typically come with desktop software or apps
to manipulate the images and for example share to facebook.

It's fun to explore doing this without relying on the official software. 

## Storage

The dokicam stores its photos and videos on its memory card in JPG and MP4
format, easily accessible via USB storage without even removing the card.

## Projection conversion

Those images and video's show the 'double fish-eye' nature of the device.
Services like Facebook, however require, require 306 imagery to be mapped
using the Equirectangular Projection. This can be achieved with `ffmpeg`
using 2 'mapping files' for your image type.

## Mapping generation

I did not find a suitable mapping for my camera online. However I did find
`projection.c` by Floris Sluiter which could generate such mapping files
for single-fisheye sources, and modified it to support double-fisheye.

Compile the generator code:

    gcc -o projection projection.c -lm

Create mapping files for video and photo's:

    ./projection -x xmap_dokicam_video.pgm -y ymap_dokicam_video.pgm -h 1440 -w 2880 -r 1440 -c 2880 -b 35 -m double
    ./projection -x xmap_dokicam.pgm -y ymap_dokicam.pgm -h 2048 -w 4096 -r 2048 -c 4096 -b 75 -m double

## Usage

Once you have created (or downloaded) the mapping files, use them with ffmpeg:

    ffmpeg -i photo.jpg -i xmap_dokicam.pgm -i ymap_dokicam.pgm -filter_complex remap out.jpg
    ffmpeg -i movie.mp4 -i xmap_dokicam_video.pgm -i ymap_dokicam_video.pgm -filter_complex remap out.mp4

For images, add exif metadata to help e.g. Facebook understand this is 360:

    exiftool -ProjectionType="equirectangular" out.jpg
