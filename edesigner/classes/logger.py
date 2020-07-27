# -*- coding: utf-8 -*-
# logger
# Jose Alfredo Martin

version = 'logger.v.9.0.0'

# Python Modules
import time
import os


class Logger:
    def __init__(self, filename, prepend_timestamp=False):
        """Initiallizes the logger instance
        filename : str (filename of the log file)
        retunrs : None"""
        self.filename = filename
        self.strtimestamp = ''
        if prepend_timestamp:
            self.strtimestamp = time.strftime("%Y%m%d", time.gmtime())
            head, tail = os.path.split(self.filename)
            self.filename = os.path.join(head, self.strtimestamp + tail)
        with open(self.filename, 'w') as f:
            f.write('log file created on: ' + time.asctime() + '\n')

    def update(self, text, to_file=True, to_screen=True):
        """updates the file
        text : str (text to print)
        to_file : bool (whether to print the text to the instance file)
        to_screen : bool (whether to print the text to the screen)
        returns : None"""
        if to_file:
            with open(self.filename, 'a') as f:
                f.write(text + '\n')
        if to_screen:
            print(text)
        return None

    def insert_file_in_log(self, filename, caption):
        """Inserts a caption and then a parameters file in the log
        filename : str (path and filename for the file)
        caption : str (Title to be placed before the file starts
        returns : None"""
        with open(filename, 'r') as f:
            lineas = f.readlines()
        self.update('', to_screen=False)
        self.update(caption, to_screen=False)
        self.update('', to_screen=False)
        for linea in lineas:
            self.update(linea.rstrip('\r\n'), to_screen=False)
        self.update('', to_screen=False)


if __name__ == '__main__':
    print(version)
