import numpy as np
import torch
import matplotlib.pyplot as plt
import matplotlib.animation as animation
plt.style.use('ggplot')

def change_colors(table, data, fgcolor, bgcolor, text_color=None, bold_rows=None, special_color_rows=None, debug=False):
    # bold_rows should be in format: {id: color}, i.e.: {0: 'blue', 4: 'yellow'}
    # special_color_rows should be in format: {id: color}, i.e.: {0: 'blue', 4: 'yellow'}
    # Changing cell colors
    colors = [['lightgray' for _ in range(len(data[0]))] for _ in range(len(data))]

    # Changing text color
    for key, cell in table.get_celld().items():
        if debug: print(f"key: {key}, cell: {cell}")
        cell.set_edgecolor(fgcolor)
        cell.set_facecolor(bgcolor)
        # case for bold_rows
        if bold_rows and (key[0] in bold_rows):
            cell.set_text_props(fontweight='bold', color=bold_rows[key[0]])
            if bold_rows: # By default bold row == Header cell
                cell.set_text_props()
                continue

        if special_color_rows and (key[0] in special_color_rows.keys()):
            if debug: print(f"special_color_rows[key[0]]: {special_color_rows[key[0]]}")
            if (special_color_rows[key[0]] is not None):
                if debug: print(f"applying {special_color_rows[key[0]]} to {cell}")
                cell.set_text_props(color=special_color_rows[key[0]])
                continue

        if text_color: # if text_color defined
            cell.set_text_props(color=text_color)
        else: # otherwise take just fgcolor
            cell.set_text_props(color=fgcolor)


class Word():
    def __init__(self, word, special_color=None, is_bold=False, positions_in={}, positions_out={}):
        """
        Word object - defines its features and positions in graph flow.

        Args:
            word (str): name of the word
            special_color (str, optional): color used in presentation of that word. Defaults to None.
            is_bold (bool, optional): weight of the word in presentation. Defaults to False.
            positions_in (dict of tuples, optional):
                        Defines positions of the Word in the flow.
                        They refer to **entrance** positions of that Word in particular tables.
                        Key should be name of table for example, while value should be tuple of x, y position of **entrance** in that table.
                        For example: {'vocab_table': (0.02, 0.65), 'embeddings_table': (0.45, 0.65)}.
                        Defaults to None.
            positions_out (dict of tuples, optional):
                        Defines positions of the Word in the flow.
                        They refer to **exit** positions of that Word in particular tables.
                        Key should be name of table for example, while value should be tuple of x, y position of **exit** in that table.
                        For example: {'vocab_table': (0.12, 0.65), 'embeddings_table': (0.55, 0.65)}.
                        Defaults to None.
        """
        self.word = word
        self.special_color = special_color
        self.is_bold = is_bold
        self.positions_in = positions_in
        self.positions_out = positions_out

    def __str__(self):
        return f"word: '{self.word}', special_color: '{self.special_color}', is_bold: '{self.is_bold}', positions_in: '{self.positions_in}', 'positions_out: '{self.positions_out}'"

    def set_special_color(self, special_color):
        self.special_color = special_color

    def set_boldness(self, is_bold):
        self.is_bold = is_bold

    def add_positions(self, name, position_in=None, position_out=None):
        # name: obj_name - name of the object the position refers to
        # position_in: (x, y) - physical positions of the word, entrance
        # position_out: (x, y) - physical positions of the word, exit
        if not position_in and not position_out:
            raise(f"At least one of 'position_in', 'position_out' must not be None. Provided  position_in: {position_in}, position_out: {position_out}.")

        if position_in:
            self.positions_in[name] = position_in
        if position_out:
            self.positions_out[name] = position_out

class Table():
    def __init__(self, title, main_pos, table_data, fontsize=7, title_fontsize=10) -> None:
        self.title = title
        self.main_pos = main_pos
        self.fontsize = fontsize
        self.title_fontsize = title_fontsize
        self.table_data = table_data
        self.colLabels = None
        self.rowLabels = None
        self.ax = None
        self.ax_table = None
        self.ax_title = None
        self.x_ratio = None
        self.y_ratio = None

        # space that table ocupies
        self.left_limit = None
        self.right_limit = None
        self.top_limit = None
        self.bottom_limit = None

    def __str__(self) -> str:
        return f"""title: {self.title},
                   main_pos: {self.main_pos},
                   x_ratio: {self.x_ratio},
                   y_ratio: {self.y_ratio},
                   fontsize: {self.fontsize},
                   title_fontsize: {self.title_fontsize},
                   colLabels: {self.colLabels},
                   rowLabels: {self.rowLabels},
                   limits (left, right, up, down): {self.left_limit}, {self.right_limit}, {self.top_limit}, {self.bottom_limit},
                   data rows number: {len(self.table_data)},
                   data columns number: {len(self.table_data[0])},
                   data: {self.table_data}"""


    def draw(self, ax, x_ratio, y_ratio, with_title=True, title_position=None, *args, **kwargs):
        x_eps = 0.05
        y_eps = 0.02
        x_ratio = x_eps*len(self.table_data[0])
        y_ratio = y_eps*len(self.table_data)
        # Calculate the positions of the borders of the table
        self.left_limit = self.main_pos[0]
        self.right_limit = self.main_pos[0] + x_ratio
        self.bottom_limit = self.main_pos[1]
        self.top_limit = self.main_pos[1] + y_ratio 

        if 'colLabels' in kwargs:
            self.colLabels = kwargs['colLabels']
        if 'rowLabels' in kwargs:
            self.rowLabels = kwargs['rowLabels']

        context_table = ax.table(self.table_data, loc='top left', cellLoc='center', bbox=[self.main_pos[0], self.main_pos[1], x_ratio, y_ratio], *args, **kwargs)
        context_table.auto_set_font_size(False)
        context_table.set_fontsize(self.fontsize)

        if with_title:
            if title_position:
                ax_title = ax.text(title_position[0], title_position[1], self.title, weight='bold', fontsize=self.title_fontsize)
            else:
                ax_title = ax.text(self.left_limit, self.top_limit+0.01, self.title, weight='bold', fontsize=self.title_fontsize)
            self.ax_title = ax_title
        # Keep object reference for future manipulations
        self.ax = ax
        self.ax_table = context_table
        self.x_ratio = x_ratio
        self.y_ratio = y_ratio


    def change_colors(self, *args, **kwargs):
        if not self.ax_table:
            raise ValueError("ax.table not defined for this Table obj, yet. Run self.draw() first.")
        if 'special_color_rows' in kwargs and self.colLabels: # move row id's by one due to header
            special_color_rows={}
            for id, color in kwargs['special_color_rows'].items():
                special_color_rows[(id+1)]=color
            kwargs['special_color_rows'] = special_color_rows

        # the same as above for bold_rows - move by one
        if 'bold_rows' in kwargs and self.colLabels: # move row id's by one due to header
            bold_rows={}
            for id, color in kwargs['bold_rows'].items():
                bold_rows[(id+1)]=color
            kwargs['bold_rows'] = bold_rows

        # table
        change_colors(self.ax_table, self.table_data, *args, **kwargs)
        # title
        if self.ax_title:
            self.ax_title.set_color(args[0])

    def row_height(self):
        # y center
        height = self.y_ratio
        if self.colLabels:
            row_height = height/(len(self.table_data)+1)
        else:
            row_height = height/len(self.table_data)
        return row_height

    def row_pos(self, row_id, pos='center'):
        # choose from pos=['center', 'left', 'right']
        # x center
        x_middle = (self.right_limit - self.left_limit)/2 + self.left_limit

        # y center
        height = self.y_ratio
        if self.colLabels:
            row_height = height/(len(self.table_data)+1)
            y_middle = self.top_limit - row_height/2 - (row_id+1)*row_height
        else:
            row_height = height/len(self.table_data)
            y_middle = self.top_limit - row_height/2 - (row_id)*row_height
        
        if pos == 'left':
            return (self.left_limit, y_middle)
        elif pos == 'right':    
            return (self.right_limit, y_middle)
        elif pos == 'center':   
            return (x_middle, y_middle)
        else:
            raise ValueError(f"pos='{pos}' not recognized. Choose from ['center', 'left', 'right']")

def generate_vocab_table_data(context_words, text_dict, embed_dims, no_filling=False, debug=False):
    torch.manual_seed(69)
    vocab_table_data = []
    embed_table_data = []
    # first two
    if no_filling is False:
        vocab_table_data.append(text_dict[0])
        vocab_table_data.append(text_dict[1])
        embed_table_data.append(torch.randn([embed_dims]).numpy())
        embed_table_data.append(torch.randn([embed_dims]).numpy())

    if debug:
        print("DEBUG (generate_vocab_table_data): text_dict: ", text_dict)
        print("DEBUG (generate_vocab_table_data): context_words len: ", len(context_words))
        print(f"[context_word.word for context_word in context_words]: {[context_word.word for context_word in context_words]}")

    #context_words_sorted_by_id = [context_word.word for context_word in context_words]
    context_words_sorted_by_id = [item for item in text_dict if item[0] in [context_word.word for context_word in context_words]]
    for context_word in context_words_sorted_by_id:
        for item in text_dict:
            if item[0] == context_word[0]:
                context_word = item
        if debug:
            print("DEBUG (generate_vocab_table_data): context word row: ", context_word)
        if context_word in vocab_table_data:
            continue

        if vocab_table_data and (int(vocab_table_data[-1][1]) != (int(context_word[1])-1)) and (no_filling is False):
            vocab_table_data.append(['...', '...'])
            embed_table_data.append(['...' for _ in range(embed_dims)])
        vocab_table_data.append(context_word)
        embed_table_data.append(torch.randn([embed_dims]).numpy())

    # last 2 rows
    if (no_filling is False) and (text_dict[-2] not in vocab_table_data):
        vocab_table_data.append(['...', '...'])
        embed_table_data.append(['...' for _ in range(embed_dims)])
        vocab_table_data.append(text_dict[-2])
        embed_table_data.append(torch.randn([embed_dims]).numpy())
    if (no_filling is False) and (text_dict[-1] not in vocab_table_data):
        vocab_table_data.append(text_dict[-1])
        embed_table_data.append(torch.randn([embed_dims]).numpy())

    if debug: print("DEBUG (generate_vocab_table_data): table data", vocab_table_data)
    return vocab_table_data, embed_table_data

def find_word_in_table(word, table, colId=0, debug=False):
    # returns id of that word in table - it's duck-typing kind of search - one should not overdepend on it - it doesn't make sense to keep it in Table class, because Table itself has no Word associations and shouldn't
    # it could be moved to Word class actually
    if type(word) is not Word:
        raise TypeError(f"Word {word} is expected to be of type Word, but found {type(word)}")
    for table_data_id, table_row in enumerate(table.table_data):
        if type(table_row) is list:
            table_row = table_row[colId]
        if debug: print(f"[DEBUG] (find_word_in_table) table_row: {table_row}, word: {word.word}")
        if table_row == word.word:
            if debug: print(f"[DEBUG] (find_word_in_table) match, returning: {table_data_id}")
            return table_data_id
    return None


def draw_pipeline(
    text,
    word_start_id,
    context_window_size=7,
    embed_dims=4,
    text_dict=None,
    ax=None,
    fgcolor='grey',
    bgcolor='black',
    base_x_position=0.02,
    base_y_position=0.5,
    between_tables_distance=0.07,
    ):

    if ax is None:
        fig, ax = plt.subplots(figsize=(16, 9))
        fig.set_facecolor(bgcolor)
        ax.axis('off')

    colors_palete = ['tab:blue', 'tab:pink', 'tab:brown', 'tab:olive', 'tab:orange', 'tab:cyan', 'tab:red', 'tab:green']

    #### Prepare data
    text_split = list(text.split())
    if not text_dict:
        text_dict = [ [word, str(key)] for key, word in enumerate(sorted(set(text_split), key=lambda v: v.upper()))]

    #sentence = list('We not are accounted poor not citizens,'.split())
    sentence = text_split[word_start_id:word_start_id+context_window_size]
    mid_position = int(len(sentence)/2)
    context_window = [Word(word) for word in sentence]
    context_words = context_window.copy()
    mid_word = context_words.pop(mid_position)

    for id, word in enumerate(context_words):
        word.set_special_color(colors_palete[id%len(colors_palete)])


    #### Draw tables and arrows

    tables_layers = {}
    # TABLE Context table
    
    # text box
    text_presentation = [[item] for item in text_split]
    row_height = 0.02
    text_field = Table("", (base_x_position, base_y_position-row_height*(len(text_split)-context_window_size-word_start_id)), text_presentation)
    text_field.draw(ax, 0.1, 0.1, edges='open')
    text_field.change_colors(fgcolor, bgcolor)

    # prepare data table/s
    context_table_data = [[item] for item in sentence]
    # draw table
    table_name = "Context\nwindow"
    context_table = Table(table_name, (base_x_position, base_y_position), context_table_data)
    context_table.draw(ax, 0.07, 0.15, title_position=(base_x_position-0.1, base_y_position+0.12))
    context_table_scr = {find_word_in_table(word, context_table, colId=0): word.special_color for word in context_words}
    context_table.change_colors(fgcolor, bgcolor, special_color_rows=context_table_scr, bold_rows={mid_position: 'grey'})
    tables_layers[table_name] = context_table

    # TABLE Vocab table
    # prepare data table/s
    vocab_words = context_words + [mid_word]
    vocab_table_data, embed_layer_data = generate_vocab_table_data(vocab_words, text_dict, embed_dims)
    # draw table
    table_name = "Vocabulary"
    vocab_table = Table(table_name, (tables_layers["Context\nwindow"].right_limit+between_tables_distance, tables_layers["Context\nwindow"].bottom_limit-0.1), vocab_table_data)
    vocab_table.draw(ax, 0.1, 0.29, colLabels=['word', 'id'])
    vocab_table_scr = {find_word_in_table(word, vocab_table, colId=0): word.special_color for word in context_words}
    vocab_table.change_colors(fgcolor, bgcolor, special_color_rows=vocab_table_scr, bold_rows={find_word_in_table(mid_word, vocab_table, colId=0): 'grey'})
    tables_layers[table_name] = vocab_table

    # ARROWS context window ranges
    row_width = 0.5
    arrow_props = dict(facecolor=mid_word.special_color, edgecolor=mid_word.special_color, arrowstyle=f'-[, widthB={row_width}')
    ax.annotate('', 
                xytext=(context_table.left_limit-0.02, context_table.row_pos(mid_position, pos='left')[1]), 
                xy=(context_table.left_limit, context_table.row_pos(mid_position, pos='left')[1]), 
                xycoords='axes fraction', 
                textcoords='axes fraction', 
                arrowprops=arrow_props)
    ax.text(context_table.left_limit-0.04, context_table.row_pos(mid_position, pos='left')[1], 'mid word',
        horizontalalignment='center',
        verticalalignment='center',
        fontsize=8, color='grey',
        transform=ax.transAxes)
    # upper context
    arrow_props = dict(facecolor=fgcolor, edgecolor=fgcolor, arrowstyle=f'-[, widthB={row_width*len(context_words)/2}')
    ax.annotate('', 
                xytext=(context_table.left_limit-0.01, context_table.top_limit - context_table.row_height()*(context_window_size-1)/4), 
                xy=(context_table.left_limit, context_table.top_limit - context_table.row_height()*(context_window_size-1)/4), 
                xycoords='axes fraction', 
                textcoords='axes fraction', 
                arrowprops=arrow_props)
    ax.text(context_table.left_limit-0.03, context_table.top_limit - context_table.row_height()*(context_window_size-1)/4, 'context',
        horizontalalignment='center',
        verticalalignment='center',
        fontsize=8, color=fgcolor,
        transform=ax.transAxes)
    # lower context
    ax.annotate('', 
                xytext=(context_table.left_limit-0.01, context_table.bottom_limit + context_table.row_height()*(context_window_size-1)/4), 
                xy=(context_table.left_limit, context_table.bottom_limit + context_table.row_height()*(context_window_size-1)/4), 
                xycoords='axes fraction', 
                textcoords='axes fraction', 
                arrowprops=arrow_props)
    ax.text(context_table.left_limit-0.03, context_table.bottom_limit + context_table.row_height()*(context_window_size-1)/4, 'context',
        horizontalalignment='center',
        verticalalignment='center',
        fontsize=8, color=fgcolor,
        transform=ax.transAxes)


    # ARROW context table -> vocab table
    vocab_words_y_positions = []
    for id, word in enumerate(context_window):
        arrow_props = dict(facecolor=word.special_color, edgecolor=word.special_color, arrowstyle='->', connectionstyle='arc, angleA=0, armA=30, angleB=0, armB=-30,rad=15')
        vocab_pos = vocab_table.row_pos(find_word_in_table(word, vocab_table, colId=0), pos='left')
        ax.annotate('', 
                    xytext=context_table.row_pos(id, pos='right'), 
                    xy=vocab_pos, 
                    xycoords='axes fraction', 
                    textcoords='axes fraction', 
                    arrowprops=arrow_props)

        # vocab table -> embed table
        if word != mid_word:
            vocab_words_y_positions.append(vocab_pos[1])

    # TABLE Embed table
    # draw table
    table_name = "Embeddings"
    embed_table = Table(table_name, (tables_layers["Vocabulary"].right_limit+between_tables_distance, tables_layers["Vocabulary"].bottom_limit), embed_layer_data)
    embed_table.draw(ax, 0.2, 0.29, colLabels=[f"dim_{id}" for id in range(len(embed_layer_data[0]))], rowLabels=[item[1] for item in vocab_table_data])
    embed_table_scr = {find_word_in_table(word, vocab_table, colId=0): word.special_color for word in context_words}
    embed_table.change_colors(fgcolor, bgcolor, special_color_rows=embed_table_scr)
    tables_layers[table_name] = embed_table


    # ARROWS vocab table -> embed table
    for id, word in enumerate(context_words):
        arrow_props = dict(facecolor=word.special_color, edgecolor=word.special_color, arrowstyle='->')
        ax.annotate('', 
                    xytext=(vocab_table.right_limit, vocab_words_y_positions[id]), 
                    xy=(embed_table.left_limit - 0.01, vocab_words_y_positions[id]), 
                    xycoords='axes fraction', 
                    textcoords='axes fraction', 
                    arrowprops=arrow_props)

    ## TABLES Embed output tables
    # prepare data table/s
    embed_output_data = []
    for word in context_words:
        data = embed_layer_data[find_word_in_table(word, vocab_table, colId=0)]
        embed_output_data.append(
            np.transpose([data])
        )
    embed_output_tables = []
    # draw tables
    table_name = "Embed Output"
    embed_out_tbl_title_fontsize = 6
    for id, embed_output_data_word in enumerate(embed_output_data):
        #y_pos = 0.84-0.1*id
        embed_table_center = (embed_table.bottom_limit + embed_table.top_limit)/2
        size_of_single_table = embed_table.row_height()*embed_dims+embed_out_tbl_title_fontsize*0.006
        top = embed_table_center + size_of_single_table*(len(embed_output_data)/2 - 1)
        y_pos = top-size_of_single_table*(id)
        embed_output_table = Table(f"\"{context_words[id].word}\"", 
                                (tables_layers["Embeddings"].right_limit+between_tables_distance, y_pos), 
                                embed_output_data_word, 
                                title_fontsize=embed_out_tbl_title_fontsize)
        embed_output_table.draw(ax, 0.1, 0.02, with_title=True)
        embed_output_table_scr = {row: context_words[id].special_color for row in range(len(embed_output_data_word))}
        embed_output_table.change_colors(fgcolor, bgcolor, special_color_rows=embed_output_table_scr)
        embed_output_tables.append(embed_output_table)
    tables_layers[table_name] = embed_output_tables
        

    # ARROWS Embeddings -> Embedding Outputs
    ## for each embed
    for id, word in enumerate(context_words):
        arrow_props = dict(facecolor=word.special_color, edgecolor=word.special_color, arrowstyle='->', connectionstyle='arc, angleA=0, armA=30, angleB=0, armB=-30, rad=15')
        y_pos_output_embed = (embed_output_tables[id].bottom_limit + embed_output_tables[id].top_limit)/2
        ax.annotate('', 
                    xytext=(embed_table.right_limit, vocab_words_y_positions[id]), 
                    xy=(embed_output_tables[id].left_limit, y_pos_output_embed), 
                    xycoords='axes fraction', 
                    textcoords='axes fraction',
                    arrowprops=arrow_props)


    # TABLE Output table
    # reuse same function, because of the shape, but rest is irrelevant
    _, output_layer_data = generate_vocab_table_data(vocab_words, text_dict, 1)
    table_name = "Output"
    output_table = Table(table_name, (tables_layers["Embed Output"][0].right_limit+2*between_tables_distance, tables_layers["Vocabulary"].bottom_limit), output_layer_data)
    output_table.draw(ax, 0.2, 0.29, colLabels=[f"prediction"], rowLabels=[item[1] for item in vocab_table_data])
    output_table_scr = {find_word_in_table(mid_word, vocab_table, colId=0): mid_word.special_color}
    output_table.change_colors(fgcolor, bgcolor, bold_rows={find_word_in_table(mid_word, vocab_table, colId=0): "grey"})
    tables_layers[table_name] = output_table

    # ARROWS Embedding Outputs -> Outputs
    ## for each embed output row to each output row
    for id, word in enumerate(context_words):
        arrow_props = dict(facecolor=word.special_color, edgecolor=word.special_color, arrowstyle='->', connectionstyle='arc, angleA=0, armA=30, angleB=0, armB=-30, rad=20')
        y_pos_output_embed = output_table.top_limit - 0.01
        for row in range(len(embed_output_tables[id].table_data)):
            for output_row in range(len(output_table.table_data)):
                x_pos, y_pos = output_table.row_pos(output_row, pos='left')
                x_pos -= 0.01
                ax.annotate('', 
                            xytext=embed_output_tables[id].row_pos(row, pos="right"), 
                            xy=(x_pos, y_pos), 
                            xycoords='axes fraction', 
                            textcoords='axes fraction',
                            arrowprops=arrow_props)

    # TABLE Expected table
    expected_layer_data = []
    for word in output_layer_data:
        if word != ['...']:
            if word == output_layer_data[find_word_in_table(mid_word, vocab_table, colId=0)]: # mid_word
                expected_layer_data.append([1.0])
            else:
                expected_layer_data.append([0.0])
        else:
            expected_layer_data.append(['...'])
    expected_layer_data = np.array(expected_layer_data)
    expected_layer_data[find_word_in_table(mid_word, vocab_table, colId=0)] = 1
    table_name = f"Expected Output:\nOne-Hot(\"{mid_word.word}\")"
    vertical_pos = 2*tables_layers["Output"].bottom_limit - tables_layers["Output"].top_limit - 0.06   
    expected_table = Table(table_name, 
                           (tables_layers["Embed Output"][0].right_limit+2*between_tables_distance, vertical_pos), 
                           expected_layer_data)
    expected_table.draw(ax, 0.2, 0.29, colLabels=[f"expected"], rowLabels=[item[1] for item in vocab_table_data])
    expected_table_scr = {find_word_in_table(mid_word, vocab_table, colId=0): mid_word.special_color}
    expected_table.change_colors(fgcolor, bgcolor, bold_rows={find_word_in_table(mid_word, vocab_table, colId=0): "grey"})
    tables_layers[table_name] = expected_table

    # calc positions for Expected and Outputs
    prediction_pos = output_table.row_pos(find_word_in_table(mid_word, vocab_table, colId=0), pos='right')
    expected_x_pos, expected_y_pos = expected_table.row_pos(find_word_in_table(mid_word, vocab_table, colId=0), pos='left')
    expected_x_pos -= 0.01
    expected_pos = (expected_x_pos, expected_y_pos)

    # ARROW Vocab(mid_word) -> Expected
    # first part
    arrow_props = dict(facecolor=mid_word.special_color, edgecolor=mid_word.special_color, arrowstyle='-|>', connectionstyle='arc, angleA=0, armA=20, angleB=0, armB=-560, rad=20')
    ax.annotate('', 
                xytext=(vocab_table.right_limit, vocab_table.row_pos(find_word_in_table(mid_word, vocab_table, colId=0)+0.02, pos='left')[1]), 
                xy=expected_pos, 
                xycoords='axes fraction', 
                textcoords='axes fraction', 
                arrowprops=arrow_props)

    # LOSS
    loss_position = (output_table.right_limit + 2*between_tables_distance, (output_table.top_limit + expected_table.bottom_limit)/2)
    ax.text(loss_position[0]+0.03, loss_position[1], 'LOSS',
        horizontalalignment='center',
        verticalalignment='center',
        fontsize=14, color='grey',
        transform=ax.transAxes)

    # ARROW Outputs -> LOSS
    # loss_position - between Expected and Outputs on the right
    loss_position = (output_table.right_limit + 2*between_tables_distance, (output_table.top_limit + expected_table.bottom_limit)/2)
    prediction_right_center_pos = (output_table.right_limit+0.01, (output_table.top_limit + output_table.bottom_limit)/2)
    expected_right_center_pos = (expected_table.right_limit+0.01, (expected_table.top_limit + expected_table.bottom_limit)/2)
    table_height_eps = 3
    arrow_props = dict(facecolor=mid_word.special_color, edgecolor=mid_word.special_color, arrowstyle=f']->, widthA={expected_table.y_ratio*(len(expected_table.table_data)+table_height_eps)}', connectionstyle='arc, angleA=0, armA=20, angleB=0, armB=-30, rad=15')
    ax.annotate('', 
                xytext=prediction_right_center_pos, 
                xy=loss_position, 
                xycoords='axes fraction', 
                textcoords='axes fraction', 
                arrowprops=arrow_props)
    ax.annotate('', 
                xytext=expected_right_center_pos, 
                xy=loss_position, 
                xycoords='axes fraction', 
                textcoords='axes fraction', 
                arrowprops=arrow_props)

    return ax


if __name__ == "__main__":
    text = '''We are accounted poor citizens, the patricians good.
    What authority surfeits on would relieve us: if they
    would yield us but the superfluity, while it were
    wholesome, we might guess they relieved us humanely;
    but they think we are too dear: the leanness that
    afflicts us, the object of our misery, is as an
    inventory to particularise their abundance; our
    sufferance is a gain to them Let us revenge this with
    our pikes, ere we become rakes: for the gods know I
    speak this in hunger for bread, not in thirst for revenge.''' 
    

    BGCOLOR="black"
    fig, ax = plt.subplots(figsize=(16, 9))

    def update(frame):
        plt.cla()
        fig.set_facecolor(BGCOLOR)
        ax.axis('off')
        draw_pipeline(text=text, 
                      word_start_id=frame, 
                      context_window_size=5, 
                      embed_dims=3, 
                      ax=ax, 
                      bgcolor=BGCOLOR)

    ani = animation.FuncAnimation(fig=fig, func=update, frames=4, interval=100)
    plt.show()