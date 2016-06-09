from os.path import join
import sys

def nc(slashcmd, val):
    return r"\newcommand{\%s}{%s}" % (slashcmd, val)

def main():
    path = sys.argv[1]

    picks = eval(sys.argv[2])

    conv = ["i", "ii", "iii", "iv", "v"]

    tex_dump = ""

    for ai in range(3):
        a_picks = picks[ai]
        a_dir = join(path, "a%d" % ai)
        ca = conv[ai]

        i = 0
        with open(join(a_dir, "conv.txt"), 'r') as f:
            for line in f:
                frame, _, F0, cov = line.split()

                if int(frame) == a_picks[i]:

                    ci = conv[i]

                    tex_dump += nc("%sfig%s" % (ca, ci), frame) + "\n"
                    tex_dump += nc("%sfz%s" % (ca, ci), F0) + "\n"
                    tex_dump += nc("%sr%s" % (ca, ci), round(float(cov), 3)) + "\n\n"

                    i += 1

                    if i == 5:
                        break

        tex_dump += "\n"

    with open("/tmp/figure_selection.tex", 'w') as f:
        f.write(tex_dump)


if __name__ == "__main__":
    main()