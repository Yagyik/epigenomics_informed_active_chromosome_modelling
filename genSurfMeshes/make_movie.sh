ffmpeg -r $3 -start_number $4 -i $1/image%04d.png -vcodec libx264 $1/$2.mp4
