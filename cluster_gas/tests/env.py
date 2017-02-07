import sys
import os

# Append module root directory to sys.path
sys.path.append(
    os.path.dirname(
        os.path.dirname(
            os.path.abspath(__file__)
        )
    )
)

# See this thread Re:tests folder. http://stackoverflow.com/a/23386287/6731049
