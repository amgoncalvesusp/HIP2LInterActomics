"""Convenience launcher: `python run.py`"""
from multiprocessing import freeze_support

freeze_support()

from luna_gui.main import main
import sys
sys.exit(main())
