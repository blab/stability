import os



os.system("python src/H3N2_process.py -v1 -y30 --stop mutations")

os.system("python feed_foldx.py")