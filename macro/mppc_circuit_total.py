#!/usr/bin/env python3
import os
import schemdraw
import schemdraw.elements as elm
schemdraw.config(font='Times New Roman')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['mathtext.rm'] = 'Times New Roman'
plt.rcParams['mathtext.it'] = 'Times New Roman:italic'
plt.rcParams['mathtext.bf'] = 'Times New Roman:bold'

# ---------- utility: rectangle using 4 lines (old schemdraw safe) ----------
def box(d, x, y, w, h, lw=1.2):
    d.add(elm.Line().at((x, y)).right().length(w).linewidth(lw))
    d.add(elm.Line().at((x+w, y)).up().length(h).linewidth(lw))
    d.add(elm.Line().at((x+w, y+h)).left().length(w).linewidth(lw))
    d.add(elm.Line().at((x, y+h)).down().length(h).linewidth(lw))

# ---------- MPPC equivalent, drawn as "rotated 90° CCW" style ----------
def draw_mppc_equiv_rot_ccw(d, at=(0, 0), w=6.0, h=10.0, n_cells=4):
    x0, y0 = at
    box(d, x0, y0, w, h, lw=1.2)

    d.add(elm.Label().at((x0 + w/2, y0 + h +1.0)).label('MPPC', loc='center'))

    # rail x positions
    left_x  = x0 + 0.22*w
    right_x = x0 + 0.78*w

    # ---- decide microcell y positions FIRST ----
    mid_y = y0 + 0.50*h
    span = 0.55*h  # microcell stack height inside the box (tweak 0.45~0.70)
    if n_cells <= 1:
        ys = [mid_y]
    else:
        start = mid_y - span/2
        step = span/(n_cells - 1)
        ys = [start + i*step for i in range(n_cells)]

    # ---- rail length matched to microcell region ----
    pad = 0.0  # extra margin above/below the first/last branch (tweak 0.2~0.6)
    rail_bot = min(ys) - pad
    rail_top = max(ys) + pad

    d.add(elm.Line().at((left_x,  rail_bot)).to((left_x,  rail_top)))
    d.add(elm.Line().at((right_x, rail_bot)).to((right_x, rail_top)))

    # bus nodes (center)
    return_node = (left_x,  mid_y)
    bias_node   = (right_x, mid_y)
    #d.add(elm.Dot().at(return_node))
    #d.add(elm.Dot().at(bias_node))

    # ---- draw microcells between rails (diode flipped + resistor on the right) ----
    lead = 0.25
    gap  = (right_x - left_x) - 2*lead
    d_len = 0.40 * gap
    r_len = 0.60 * gap

    for yi in ys:
        d.add(elm.Line().at((left_x, yi)).right().length(lead))

        di = elm.Diode().right().length(d_len).reverse()
        d.add(di)

        d.add(elm.Resistor().right().length(r_len))
        d.add(elm.Line().right().length(lead))  # ends right at the rail

    d.add(elm.Label().at((x0 + w/2, y0 + h/2)).label('...', loc='center'))

    # signal out to the right (clean connect to PZ)
    signal_out = (x0 + w + 1.0, mid_y)
    d.add(elm.Line().at(bias_node).right().to(signal_out))

    return bias_node, return_node, signal_out

# ---------- PZ cancellation network ----------
def draw_pz_network(d, in_pt, Rser=r'$22\,\Omega$', Cser='100 pF', Rpar=r'$390\,\Omega$',
                    dx=7.0, dy=2.2):
    # lead from input node
    d.add(elm.Line().at(in_pt).right().length(0.8))
    left = d.here
    #d.add(elm.Dot().at(left))

    # top branch: R + C
    d.add(elm.Line().at(left).up().length(dy/2))
    d.add(elm.Resistor().right().length(dx/2).label(Rser))
    d.add(elm.Capacitor().right().length(dx/2).label(Cser))
    top_end = d.here
    d.add(elm.Line().down().length(dy/2))

    # bottom branch: Rpar
    d.add(elm.Line().at(left).down().length(dy/2))
    d.add(elm.Resistor().right().length(dx).label(Rpar))
    bot_end = d.here
    d.add(elm.Line().up().length(dy/2))

    out_pt = (top_end[0], left[1])
    #d.add(elm.Dot().at(out_pt))

    # box around PZ (compact)
    bx = left[0] - 1.0
    by = left[1] - (dy/2 + 0.8)
    bw = dx + 2.0
    bh = dy + 2.0
    #box(d, bx, by, bw, bh, lw=1.0)
    d.add(elm.Label().at((bx + bw/2, by + bh + 0.5)).label('Pole-zero cancellation', loc='center'))

    return out_pt

# ---------- AD8000 inverting stage (your reference style) ----------
def draw_ad8000_stage(d, in_pt, RF_label=r'$22\,\Omega$', Cout_label='C53\n0.1u',
                      op_label='U2\nAD8000'):
    # lead into summing node
    d.add(elm.Line().at(in_pt).right().length(1.8))
    sum_node = d.here
    #d.add(elm.Dot().at(sum_node))

    op = d.add(elm.Opamp().at((sum_node[0] + 6.0, sum_node[1])).right())

    # --- clean orthogonal routing to opamp (-) input ---
    # 1) horizontal from summing node
    d.add(elm.Line().at(sum_node).right().length(1.2))
    mid = d.here

    # 2) vertical to opamp input height
    d.add(elm.Line().at(mid).up().to((mid[0], op.in2[1])))

    # 3) horizontal into opamp (-) input
    d.add(elm.Line().right().to(op.in1))
    #d.add(elm.Dot().at(op.in2))

    # (+) to ground
    d.add(elm.Line().at(op.in2).left().length(1.0))
    d.add(elm.Ground())

    # feedback loop above
    gap = 0.9          # (-)핀에서 왼쪽으로 뺄 거리
    loop_h = 2.0       # 루프 높이
    loop_w = 4.2       # 루프 폭(저항 길이)
    
    # 1) from (-) input go LEFT a bit (creates the visual gap)
    d.add(elm.Line().at(op.in1).left().length(gap))
    fb_left = d.here

    # 2) go UP
    d.add(elm.Line().at(fb_left).up().length(loop_h))
    fb_top_left = d.here

    # 3) resistor across the top
    d.add(elm.Resistor().at(fb_top_left).right().length(loop_w).label(RF_label, loc='top'))
    fb_top_right = d.here

    # 4) go DOWN to the output y (not necessarily op.out)
    out_y = op.center[1]   # output at symbol center height
    d.add(elm.Line().at(fb_top_right).down().to((fb_top_right[0], out_y)))

    # output + coupling cap
    d.add(elm.Line().at(op.out).right().length(1.2))
    d.add(elm.Capacitor().right().length(2.0).label(Cout_label, loc='top'))
    x, y = d.here
    d.add(elm.Label().at((x + 1.5, y)).label('Output', loc='center'))

    d.add(elm.Label().at((op.center[0], op.center[1] - 2.2)).label(op_label, loc='center'))
    return d.here

# ===================== main =====================
if __name__ == '__main__':
    VBIAS = 'Vbias'
    R_SER = '$22\,\Omega$'
    C_SER = '100 pF'
    R_PAR = '$390\,\Omega$'
    RF    = '$10\,\mathrm{k}\Omega$'
    COUT  = '$0.1\,\mu\mathrm{F}$'
    OPAMP = 'Op. Amp'

    d = schemdraw.Drawing()

    d.config(
        margin=0.5,
        unit=1.4,
        fontsize=25
    )

    # ---- MPPC rotated block (left) ----
    bias_node, return_node, sig_out = draw_mppc_equiv_rot_ccw(d, at=(0, 0), w=6.0, h=10.0, n_cells=4)

    # Vbias source -> return_node (instead of GND)
    vb_src2 = (return_node[0] - 4.5, return_node[1] + 3.0)   # 위치는 취향대로 조절
    d.add(elm.SourceV().at(vb_src2))
    d.add(elm.Line().at(vb_src2).down().to((vb_src2[0], return_node[1])))
    
    #d.add(elm.Dot().at(return_node))

    trimmer = (return_node[0] - 4.0, return_node[1])
    pot = d.add(elm.Potentiometer().at(trimmer).right().label('VR\n$100\\,\\mathrm{k\\Omega}$', loc='top'))

    vb_src_end = (vb_src2[0],return_node[1])
    vr1_in = pot.start
    vr1_out = pot.end

    d.add(elm.Line().at(vb_src_end).to(vr1_in))
    d.add(elm.Line().at(vr1_out).to(return_node))



    # ---- 2 MΩ after trimmer ----
    trim_resistor = (return_node[0] - 2.0, return_node[1]-1.0)
    r68 = d.add(elm.Resistor().at(trim_resistor).down().label(r'$2\,\mathrm{M\Omega}$', loc='top'))
    
    r68_top = r68.start
    r68_bot = r68.end

    trim_resistor_start = (trim_resistor[0], trim_resistor[1] +1.0)
    d.add(elm.Line().at(trim_resistor_start).to(r68_top))
    d.add(elm.Line().at(r68_bot).down().length(0.2))
    d.add(elm.Ground())
    


    d.add(elm.Line().at((vb_src2[0],vb_src2[1]+1.3)).up().length(1.0))
    d.add(elm.Ground().up())

    d.add(elm.Label().at((vb_src2[0] - 1.5, vb_src2[1])).label(r'$\mathrm{V}_{\mathrm{bias}}$', loc='right'))
    

    # ---- PZ network directly from MPPC output (already outside box) ----
    pz_out = draw_pz_network(d, sig_out, Rser=R_SER, Cser=C_SER, Rpar=R_PAR, dx=7.0, dy=2.2)

    d.add(elm.Label().at((pz_out[0] + 3.5, pz_out[1]-0.5)).label('×16 MPPCs', loc='left'))
    
    # ---- AD8000 stage ----
    draw_ad8000_stage(d, pz_out, RF_label=RF, Cout_label=COUT, op_label=OPAMP)

    

    # save
    d.save('mppc_readout.pdf')
    d.draw(show=True)
    
