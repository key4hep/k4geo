import sys

# Define variables for the number of lines before and after the match
lines_before = 3  # Number of lines to print before the match
lines_after = 2  # Number of lines to print after the match

buffer = []  # String buffer, size between 0 and 'lines_before'
skip = 0

for line in sys.stdin:
    # if pattern is found, clear buffer and set the number of lines to skip in the following if
    if "some facets have wrong orientation!" in line:
        buffer.clear()  # Clear the buffer before the matched line
        skip = lines_after  # skip the specified number of lines after the match
        continue

    # skip a number of lines after finding a match
    if skip > 0:
        skip -= 1
        continue

    buffer.append(line)

    # if buffer has enough lines, pop first and print it
    if len(buffer) > lines_before:
        sys.stdout.write(buffer.pop(0))  # Print the earliest line in the buffer
