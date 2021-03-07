from sisyphus import *
from scipy.constants import h, c


def parse_transitions(transitions):
    output = []

    for t in transitions:
        pair = []
        for state in t.split(' -> '):
            store = [
                int(state[0]),
                convert_orbital_letter_to_number(state[1]),
                float(state[2:].split('|')[0]),
                float(state[2:].split('|')[1].split(',')[0]),
                float(state[2:].split('|')[1].split(',')[1].split('>')[0])
            ]
            pair.append(store)
        output.append(pair)
    return output

def determine_index(state):
    n, l, j, f, m = state
    output = 0
    if f == j + 0.5:
        output = int(f - m)
    elif f == j - 0.5:
        output = int((2*j+2) + (f - m))
    print(output)
    return output


def get_linewidth(states):
    linewidth = 0
    label1 = ''.join((str(states[0][0]), convert_orbital_number_to_letter(states[0][1]), str(states[0][2])))
    label2 = ''.join((str(states[1][0]), convert_orbital_number_to_letter(states[1][1]), str(states[1][2])))
    query_string = '{} -> {}'.format(label1, label2)
    if query_string in nist_decay_rates.keys():
        linewidth = nist_decay_rates[query_string]
    print('linewidth', query_string, states, linewidth)
    return linewidth


def generate_resonance_plot(transitions, B, energy_offset=0, energy_scaling=1):

    parsed_transitions = parse_transitions(transitions)

    fig = plt.figure(figsize=(5, 6))
    ax = fig.add_subplot(111)
    fig.subplots_adjust(bottom=0.2)

    n = np.arange(0, 0.1, 0.0001, dtype=float)

    atoms = {}
    lines = []
    for t in parsed_transitions:
        temp_lines = []
        for s in t:
            label = ''.join([str(s[0]), str(s[1])])
            if label not in atoms.keys():
                atoms[label] = Atom(s[0], s[1], B_field=B, energy_offset=energy_offset, energy_scaling=energy_scaling)
            temp_lines.append(atoms[label].eigen_range(n)[:,determine_index(s)]/(-1*h))
            print(s, determine_index(s))
        lines.append([temp_lines[1] - temp_lines[0], get_linewidth(t)])
    
    print(len(lines))
    for line, i in zip(lines, range(len(lines))):
        plt.plot(n, line[0], label=transitions[i])
        plt.fill_between(n, line[0]-line[1]/2, line[0]+line[1]/2, alpha=0.5)

    for laser in lasers:
        line = np.full(n.shape, c/laser[0])
        plt.plot(n, line)
        plt.fill_between(n, line - laser[1]/2, line + laser[1]/2, alpha=0.5)
    plt.legend()
    plt.show()





if __name__ == '__main__':
    # Define B field 
    profile = np.array([lambda x: 0,
                        lambda y: 0,
                        lambda z : z ])
    B = Field(profile)

    # Specify transitions to plot
    transitions = [
        '2p1.5|2,2> -> 1s0.5|1,1>',
        '2p1.5|2,1> -> 1s0.5|1,0>',
        '2p1.5|2,0> -> 1s0.5|1,-1>',
    ]

    # Specify lasers to plot
    lasers = [
        # [c/(71.5e9-2.4660719e15), 0]
    ]

    generate_resonance_plot(transitions, B)

    





