
module fir_filter (
    input signed [15:0] x,      // Input value
    output signed [44:0] h0,
    output signed [44:0] h1,
    output signed [44:0] h2,
    output signed [44:0] h3,
    output signed [44:0] h4,
    output signed [44:0] h5,
    output signed [44:0] h6,
    output signed [44:0] h7,
    output signed [44:0] h8,
    output signed [44:0] h9,
    output signed [44:0] h10,
    output signed [44:0] h11,
    output signed [44:0] h12,
    output signed [44:0] h13,
    output signed [44:0] h14,
    output signed [44:0] h15,
    output signed [44:0] h16,
    output signed [44:0] h17,
    output signed [44:0] h18,
    output signed [44:0] h19,
    output signed [44:0] h20,
    output signed [44:0] h21,
    output signed [44:0] h22,
    output signed [44:0] h23,
    output signed [44:0] h24,
    output signed [44:0] h25,
    output signed [44:0] h26,
    output signed [44:0] h27,
    output signed [44:0] h28,
    output signed [44:0] h29,
    output signed [44:0] h30,
    output signed [44:0] h31
);

    // Create shifted versions of input
    wire signed [44:0] x_shift17 = x <<< 17;
    wire signed [44:0] x_shift15 = x <<< 15;
    wire signed [44:0] x_shift12 = x <<< 12;
    wire signed [44:0] x_shift11 = x <<< 11;
    wire signed [44:0] x_shift10 = x <<< 10;
    wire signed [44:0] x_shift9 = x <<< 9;
    wire signed [44:0] x_shift8 = x <<< 8;
    wire signed [44:0] x_shift7 = x <<< 7;
    wire signed [44:0] x_shift6 = x <<< 6;
    wire signed [44:0] x_shift5 = x <<< 5;
    wire signed [44:0] x_shift4 = x <<< 4;
    wire signed [44:0] x_shift3 = x <<< 3;
    wire signed [44:0] x_shift2 = x <<< 2;
    wire signed [44:0] x_shift0 = x <<< 0;

    // Cross-CSE: shared pattern "0+0+"
    wire signed [44:0] _cse_0 = x_shift8 + x_shift6;

    // h0: 00000000000000000000000+0+0-0+
    wire signed [44:0] h0 = (_cse_0 >>> 2) + -x_shift2 + x_shift0;

    // h1: 0000000000000000000+0000-00+0+
    wire signed [44:0] h1 = x_shift10 - x_shift5 + (_cse_0 >>> 6);

    // h2: 00000000000000000000000+0-0+0-
    wire signed [44:0] h2 = x_shift6 - x_shift4 + x_shift2 - x_shift0;

    // h3: 00000000000000000000+00+0+000-
    wire signed [44:0] h3 = x_shift9 + (_cse_0 >>> 2) + -x_shift0;

    // h4: 0000000000000000000000+0-00+0-
    wire signed [44:0] h4 = x_shift7 - x_shift5 + x_shift2 - x_shift0;

    // h5: 00000000000000000000+000-00+0+
    wire signed [44:0] h5 = x_shift9 - x_shift5 + (_cse_0 >>> 6);

    // h6: 00000000000000000000+000+0+00-
    wire signed [44:0] h6 = x_shift9 + (_cse_0 >>> 3) + -x_shift0;

    // h7: 0000000000000000000000+00+0+0-
    wire signed [44:0] h7 = x_shift7 + (_cse_0 >>> 4) + -x_shift0;

    // h8: 000000000000000000000+00+00+0+
    wire signed [44:0] h8 = x_shift8 + x_shift5 + (_cse_0 >>> 6);

    // h9: 00000000000000000000+000+0-00-
    wire signed [44:0] h9 = x_shift9 + x_shift5 - x_shift3 - x_shift0;

    // h10: 000000000000000000000+00-00-0+
    wire signed [44:0] h10 = x_shift8 - x_shift5 - x_shift2 + x_shift0;

    // h11: 000000000000000000000+0+00+00-
    wire signed [44:0] h11 = _cse_0 + x_shift3 - x_shift0;

    // h12: 000000000000000000+0-0000-000-
    wire signed [44:0] h12 = x_shift11 - x_shift9 - x_shift4 - x_shift0;

    // h13: 00000000000000000000+00-000-0+
    wire signed [44:0] h13 = x_shift9 - x_shift6 - x_shift2 + x_shift0;

    // h14: 0000000000000000000-000+00+00+
    wire signed [44:0] h14 = -x_shift10 + x_shift6 + x_shift3 + x_shift0;

    // h15: 00000000000000-00-000000000-0+
    wire signed [44:0] h15 = -x_shift15 - x_shift12 - x_shift2 + x_shift0;

    // h16: 000000000000000000-0+000000-0+
    wire signed [44:0] h16 = -x_shift11 + x_shift9 - x_shift2 + x_shift0;

    // h17: 00000000000000000000-0+00-000-
    wire signed [44:0] h17 = -x_shift9 + x_shift7 - x_shift4 - x_shift0;

    // h18: 0000000000000000000000-0+0+00+
    wire signed [44:0] h18 = -x_shift7 + (_cse_0 >>> 3) + x_shift0;

    // h19: 0000000000000000000-0000+00+0+
    wire signed [44:0] h19 = -x_shift10 + x_shift5 + (_cse_0 >>> 6);

    // h20: 00000000000000000000-000+0+00-
    wire signed [44:0] h20 = -x_shift9 + (_cse_0 >>> 3) + -x_shift0;

    // h21: 0000000000000000000+000-0-000-
    wire signed [44:0] h21 = x_shift10 - x_shift6 - x_shift4 - x_shift0;

    // h22: 000000000000000000000+00-00-0+
    wire signed [44:0] h22 = x_shift8 - x_shift5 - x_shift2 + x_shift0;

    // h23: 000000000000000000+0-0000-000+
    wire signed [44:0] h23 = x_shift11 - x_shift9 - x_shift4 + x_shift0;

    // h24: 00000000000000000000+00-000+0-
    wire signed [44:0] h24 = x_shift9 - x_shift6 + x_shift2 - x_shift0;

    // h25: 0000000000000000000+00-0000-0+
    wire signed [44:0] h25 = x_shift10 - x_shift7 - x_shift2 + x_shift0;

    // h26: 000000000000+0-0000000+000000-
    wire signed [44:0] h26 = x_shift17 - x_shift15 + x_shift7 - x_shift0;

    // h27: 0000000000000000000000+0+0-00+
    wire signed [44:0] h27 = (_cse_0 >>> 1) + -x_shift3 + x_shift0;

    // h28: 0000000000000000000000+00-0-0+
    wire signed [44:0] h28 = x_shift7 - x_shift4 - x_shift2 + x_shift0;

    // h29: 000000000000000000+000+0000+0+
    wire signed [44:0] h29 = x_shift11 + x_shift7 + (_cse_0 >>> 6);

    // h30: 00000000000000000000+0-0-0000+
    wire signed [44:0] h30 = x_shift9 - x_shift7 - x_shift5 + x_shift0;

    // h31: 000000000000+0000000000-0-000-
    wire signed [44:0] h31 = x_shift17 - x_shift6 - x_shift4 - x_shift0;
endmodule

