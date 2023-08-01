
import os, sys
import gzip
import shutil

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
  
SED_DIR = '../SED'

def main():
    for i, (subdir, dirs, files) in enumerate(os.walk(SED_DIR)):
        if i>0: #ignore first parent directory
            # print(i)
            print(subdir)
            print("")
            for file in files:
                # print(os.path.join(subdir, file))
                try:
                    file_name = os.path.join(subdir, file)
                    # print(file_name)
                    if file_name.endswith(".gz"):
                        print(file_name)
                        output_name = file_name[0:-3] #removed extension
                        print(output_name)
                        with gzip.open(file_name, 'rb') as f_in:
                            with open(output_name, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                        os.remove(file_name)
                except FileNotFoundError:
                    print("Trying to unzip SED files in directory SED, found none.")
                    sys.exit(1)

if __name__ == "__main__":
    main()