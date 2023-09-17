import os
import pandas as pd
import queue

class FilesInPath:

    def __initQueue(self):
        for root, dirs, files in os.walk(self.path):
            for file in files:
                try:
                    if file.endswith(".csv"):
                        rel_dir = os.path.relpath(root, self.path)
                        rel_file = os.path.join(rel_dir,file)
                        self.queue.put(rel_file)
                except:
                    continue

    def __init__(self, path):
        self.path = path
        self.queue = queue.Queue()
        self.__initQueue()

    def read_file(self, file_name):
        """
        This function is reading a parquet file contains several tweets
        The file location is given as a string as an input to this function.
        :param file_name: string - indicates the path to the file we wish to read.
        :return: a dataframe contains tweets.
        """

        df = pd.read_csv(self.path + '/' + file_name, sep=',', engine='python')
        return df

    def __iter__(self):
        return self

    def __next__(self):
        if self.queue.empty():
            raise StopIteration
        file_path = self.queue.get()
        return self.read_file(file_path), file_path

    def getNext(self):
        if self.queue.empty():
            return None
        else:
            file_path = self.queue.get()
            return self.read_file(file_path)

