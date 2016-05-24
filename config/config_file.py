#!/usr/bin/env python

# config_file.py

# By Nathan C. Hearn
#    May 2, 2007
#
# Read through ".ini"-style configuration files


# Constants

comment_marker = '#'
continuation_marker = "\\"

makefile_assign_marker = "="

help_string_marker = ":"

left_bracket_marker = "["
right_bracket_marker = "]"


class ContextSettings :

    def __init__(self, context_name, help_string=None) :

        self.__context_name = context_name

        self.__help_str = help_string

        self.__key_values = dict()


    def get_context_name(self) :
        return self.__context_name


    def get_help_string(self) :
        return self.__help_str


    def has_key(self, key) :
        return self.__key_values.has_key(key)
    

    def add_value(self, key, value) :

        # WARN USER ABOUT DUPLICATE KEYS???

        self.__key_values[key] = value


    def append(self, other_context) :

        key_list = other_context.get_keys()

        for key in key_list :

            self.add_value(key, other_context[key])


    def get_value(self, key) :

        ret_val = None

        if (self.__key_values.has_key(key)) :
            ret_val = self.__key_values[key]

        return ret_val


    def get_value_str(self, key) :

        ret_str = str()

        val = self.get_value(key)

        if (isinstance(val, list)) :
            
            num_elem = len(val)

            if (num_elem > 0) :
                
                ret_str = str(val[0])

                for idx in range(1, num_elem) :

                    ret_str += (" " + str(val[idx]))

        else :
            ret_str = str(val)

        return ret_str


    def get_keys(self) :
        return self.__key_values.keys()


    def __getitem__(self, key) :
        return self.get_value(key)


    def __setitem__(self, key, value) :
        self.add_value(key, value)


    def __iadd__(self, other_context) :
        self.append(other_context)
        

    def __str__(self) :

        if (self.__help_str is None) :
            ret_str = left_bracket_marker + " " + self.__context_name + " " \
                      + right_bracket_marker + "\n"
        else :
            ret_str = left_bracket_marker + " " + self.__context_name + " " \
                      + help_string_marker + " " + self.__help_str + " " \
                      + right_bracket_marker + "\n"

        for key in self.__key_values :

            ret_str += (key + " " + self.get_value_str(key) + "\n")

        return ret_str


# Function definitions

def remove_comments(line) :

    # Returns line with out characters after comment marker
    # ... or spaces just before marker

    work_str = line.strip()

    try :
        work_str = work_str.partition(comment_marker)[0]
    except AttributeError :
        elements = work_str.split(comment_marker)

        if (len(elements) > 0) :
            work_str = elements[0]
        else :
            work_str = ""

    return work_str.strip()


def read_lines(filename) :

    # Read the file, remove comments, and bridge continuation lines

    f = file(filename, "r")

    filelines = f.readlines()

    num_lines = len(filelines)

    current_line = str("")

    line_list = list()

    line_index = int(0)

    while (line_index < num_lines) :

        line = str("")

        still_continuing = True

        while (still_continuing) :

            if (line_index < num_lines) :
            
                current_line = remove_comments(filelines[line_index])
                line_index += 1

                num_chars = len(current_line)

                if (num_chars > 0) :

                    last_char_index = num_chars - 1

                    if (current_line[last_char_index] == continuation_marker) :

                        # Remove the marker, but not spaces before the marker

                        current_line = current_line[0:last_char_index]

                    else :

                        still_continuing = False

                    line += current_line

                else :
                    still_continuing = False

        line_list.append(line)

    return line_list


def read_config(filename) :

    # Read the file

    filelines = read_lines(filename)

    num_lines = len(filelines)

    # Build the context lists

    context_list = dict()

    current_context = ContextSettings("")

    for line in filelines :

        # Clean the line and remove any comments

        work_line = remove_comments(line)

        num_char = len(work_line)

        if (num_char > 0) :

            # Check for context switch; otherwise read key-value pair

            if (work_line[0] == left_bracket_marker) :

                # Build string from characters up to trailing right-bracket

                rbracket = work_line.find(right_bracket_marker)

                bracket_text = work_line[1:rbracket].strip()

                # If a ":" separator exists, title is before, help string after

                new_title = None
                help_str = None

                separator = bracket_text.find(help_string_marker)

                if (separator >= 0) :

                    text_len = len(bracket_text)
                    
                    new_title = bracket_text[0:separator].strip()
                    help_str = bracket_text[(separator + 1):text_len].strip()

                else :
                    
                    new_title = bracket_text

                # Move the current context the context_list and create new one

                old_title = current_context.get_context_name()

                if (context_list.has_key(old_title)) :
                    context_list[old_title].append(current_context)
                else :
                    context_list[old_title] = current_context

                current_context = ContextSettings(new_title, help_str)

            else :

                words = work_line.split()

                key_word = words[0]

                values = list()

                num_words = len(words)

                if (num_words < 2) :
                    values.append(str(""))
                else :
                    for idx in range(1, num_words) :
                        values.append(words[idx])

                current_context[key_word] = values

        # Save the last context created

        old_title = current_context.get_context_name()

        if (context_list.has_key(old_title)) :
            context_list[old_title].append(current_context)
        else :
            context_list[old_title] = current_context
        
    return context_list


def read_makefile(filename) :

    # Read the file

    filelines = read_lines(filename)

    num_lines = len(filelines)

    context = ContextSettings("", "Read from file " + filename)

    for line in filelines :

        num_char = len(line)

        if (num_char > 0) :

            keyword = None
            def_words = None

            try :
                
                keyword, sep, definition \
                         = line.partition(makefile_assign_marker)

                keyword = keyword.strip()

                def_words = definition.split()
                
            except AttributeError :
                
                line_elements = line.split(makefile_assign_marker, 1)

                keyword = line_elements[0].strip()

                num_elements = len(line_elements)

                if (num_elements > 1) :
                    def_words = line_elements[1].split()
                else :
                    def_words = list()

            num_def_words = len(def_words)

            def_list = list()

            if (num_def_words < 1) :
                def_list.append(str(""))
            else :
                def_list = def_words

            context[keyword] = def_list

    return context


if (__name__ == "__main__") :
    
    from sys import argv, stderr, stdout

    argc = len(argv)

    if (argc < 2) :
        stderr.write("Usage: " + argv[0] + " filename\n")
        exit(1)

    argPtr = int(1)

    filename = str(argv[argPtr])
    argPtr += 1

    context_list = read_config(filename)

    for key in context_list :
        print context_list[key]

    
