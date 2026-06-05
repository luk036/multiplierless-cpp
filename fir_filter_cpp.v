
module fir_filter (
    input signed [15:0] x,      // Input value
    output signed [41:0] h0
    output signed [41:0] h1
    output signed [41:0] h2
    output signed [41:0] h3
    output signed [41:0] h4
    output signed [41:0] h5
    output signed [41:0] h6
    output signed [41:0] h7
    output signed [41:0] h8
    output signed [41:0] h9
    output signed [41:0] h10
    output signed [41:0] h11
    output signed [41:0] h12
    output signed [41:0] h13
    output signed [41:0] h14
    output signed [41:0] h15
    output signed [41:0] h16
    output signed [41:0] h17
    output signed [41:0] h18
    output signed [41:0] h19
    output signed [41:0] h20
    output signed [41:0] h21
    output signed [41:0] h22
    output signed [41:0] h23
    output signed [41:0] h24
    output signed [41:0] h25
    output signed [41:0] h26
    output signed [41:0] h27
    output signed [41:0] h28
    output signed [41:0] h29
    output signed [41:0] h30
    output signed [41:0] h31
);

    // Create shifted versions of input
    wire signed [41:0] x_shift17 = x <<< 17;
    wire signed [41:0] x_shift16 = x <<< 16;
    wire signed [41:0] x_shift14 = x <<< 14;
    wire signed [41:0] x_shift13 = x <<< 13;
    wire signed [41:0] x_shift12 = x <<< 12;
    wire signed [41:0] x_shift11 = x <<< 11;
    wire signed [41:0] x_shift10 = x <<< 10;
    wire signed [41:0] x_shift9 = x <<< 9;
    wire signed [41:0] x_shift8 = x <<< 8;
    wire signed [41:0] x_shift7 = x <<< 7;
    wire signed [41:0] x_shift6 = x <<< 6;
    wire signed [41:0] x_shift5 = x <<< 5;
    wire signed [41:0] x_shift4 = x <<< 4;
    wire signed [41:0] x_shift3 = x <<< 3;
    wire signed [41:0] x_shift2 = x <<< 2;
    wire signed [41:0] x_shift0 = x <<< 0;

    // Cross-CSE: shared pattern "+0-"
    wire signed [41:0] _cse_0 = x_shift14 - x_shift12;

    // h0: 000000000000000000+0-0+000+
    wire signed [41:0] h0 = (_cse_0 >>> 6) + x_shift4 + x_shift0;

    // h1: 000000000000000000+0-0-000-
    wire signed [41:0] h1 = (_cse_0 >>> 6) + -x_shift4 - x_shift0;

    // h2: 000000000000000000+0+0-000+
    wire signed [41:0] h2 = x_shift8 + (_cse_0 >>> 8) + x_shift0;

    // h3: 0000000000000000000+00-0+0+
    wire signed [41:0] h3 = x_shift7 - x_shift4 + x_shift2 + x_shift0;

    // h4: 0000000000000000000+0+00+0-
    wire signed [41:0] h4 = x_shift7 + x_shift5 + (_cse_0 >>> 12);

    // h5: 000000000000000000+0-0+000+
    wire signed [41:0] h5 = (_cse_0 >>> 6) + x_shift4 + x_shift0;

    // h6: 00000000000000+0000-000+00+
    wire signed [41:0] h6 = x_shift12 - x_shift7 + x_shift3 + x_shift0;

    // h7: 00000000000000000+000+0+00+
    wire signed [41:0] h7 = x_shift9 + x_shift5 + x_shift3 + x_shift0;

    // h8: 0000000000000+00+000000-00+
    wire signed [41:0] h8 = x_shift13 + x_shift10 - x_shift3 + x_shift0;

    // h9: 0000000000000000000+00+0-0-
    wire signed [41:0] h9 = x_shift7 + (_cse_0 >>> 10) + -x_shift0;

    // h10: 00000000000000000+0000-0-0+
    wire signed [41:0] h10 = x_shift9 - x_shift4 - x_shift2 + x_shift0;

    // h11: 00000000000000000+0-000+00+
    wire signed [41:0] h11 = (_cse_0 >>> 5) + x_shift3 + x_shift0;

    // h12: 000000000000000000+000+0-0-
    wire signed [41:0] h12 = x_shift8 + (_cse_0 >>> 10) + -x_shift0;

    // h13: 0000000000000000+00000+0+0-
    wire signed [41:0] h13 = x_shift10 + x_shift4 + (_cse_0 >>> 12);

    // h14: 000000000000+0-0000000+000-
    wire signed [41:0] h14 = _cse_0 + x_shift4 - x_shift0;

    // h15: 0000000000000000-00+0+0000-
    wire signed [41:0] h15 = -x_shift10 + x_shift7 + x_shift5 - x_shift0;

    // h16: 000000000000000-0+00+00000-
    wire signed [41:0] h16 = -x_shift11 + x_shift9 + x_shift6 - x_shift0;

    // h17: 0000000000000000-00+00+000+
    wire signed [41:0] h17 = -x_shift10 + x_shift7 + x_shift4 + x_shift0;

    // h18: 00000000000000-00+0+000000+
    wire signed [41:0] h18 = -x_shift12 + x_shift9 + x_shift7 + x_shift0;

    // h19: 00000000000000000000-0+0+0-
    wire signed [41:0] h19 = -x_shift6 + x_shift4 + (_cse_0 >>> 12);

    // h20: 000000000000000000-000+0-0+
    wire signed [41:0] h20 = -x_shift8 + (_cse_0 >>> 10) + x_shift0;

    // h21: 00000000000000000000-0+0-0-
    wire signed [41:0] h21 = -x_shift6 + (_cse_0 >>> 10) + -x_shift0;

    // h22: 000000000000000000+00-00-0+
    wire signed [41:0] h22 = x_shift8 - x_shift5 - x_shift2 + x_shift0;

    // h23: 0000000000000000+00+0-0000+
    wire signed [41:0] h23 = x_shift10 + (_cse_0 >>> 7) + x_shift0;

    // h24: 00000000000000000000+0-0+0+
    wire signed [41:0] h24 = (_cse_0 >>> 8) + x_shift2 + x_shift0;

    // h25: 00000000000000+0000-0-0000-
    wire signed [41:0] h25 = x_shift12 - x_shift7 - x_shift5 - x_shift0;

    // h26: 000000000+000-000-00000000-
    wire signed [41:0] h26 = x_shift17 - x_shift13 - x_shift9 - x_shift0;

    // h27: 000000000000000000+0-0+000+
    wire signed [41:0] h27 = (_cse_0 >>> 6) + x_shift4 + x_shift0;

    // h28: 000000000000000000+0+000+0+
    wire signed [41:0] h28 = x_shift8 + x_shift6 + x_shift2 + x_shift0;

    // h29: 000000000000000+00-0000+00-
    wire signed [41:0] h29 = x_shift11 - x_shift8 + x_shift3 - x_shift0;

    // h30: 0000000000+000000+000-0000+
    wire signed [41:0] h30 = x_shift16 + x_shift9 - x_shift5 + x_shift0;

    // h31: 00000000000000000+000-00+0+
    wire signed [41:0] h31 = x_shift9 - x_shift5 + x_shift2 + x_shift0;
endmodule
