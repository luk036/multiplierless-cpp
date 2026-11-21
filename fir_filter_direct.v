
module fir_filter_direct (
    input signed [15:0] x,      // Input value
    output signed [40:0] h0,
    output signed [40:0] h1,
    output signed [40:0] h2,
    output signed [40:0] h3,
    output signed [40:0] h4,
    output signed [40:0] h5,
    output signed [40:0] h6,
    output signed [40:0] h7,
    output signed [40:0] h8,
    output signed [40:0] h9,
    output signed [40:0] h10,
    output signed [40:0] h11,
    output signed [40:0] h12,
    output signed [40:0] h13,
    output signed [40:0] h14,
    output signed [40:0] h15,
    output signed [40:0] h16,
    output signed [40:0] h17,
    output signed [40:0] h18,
    output signed [40:0] h19,
    output signed [40:0] h20,
    output signed [40:0] h21,
    output signed [40:0] h22,
    output signed [40:0] h23,
    output signed [40:0] h24,
    output signed [40:0] h25,
    output signed [40:0] h26,
    output signed [40:0] h27,
    output signed [40:0] h28,
    output signed [40:0] h29,
    output signed [40:0] h30,
    output signed [40:0] h31
);

    // Create shifted versions of input
    wire signed [40:0] x_shift22 = x <<< 22;
    wire signed [40:0] x_shift21 = x <<< 21;
    wire signed [40:0] x_shift20 = x <<< 20;
    wire signed [40:0] x_shift19 = x <<< 19;
    wire signed [40:0] x_shift18 = x <<< 18;
    wire signed [40:0] x_shift17 = x <<< 17;
    wire signed [40:0] x_shift16 = x <<< 16;
    wire signed [40:0] x_shift15 = x <<< 15;
    wire signed [40:0] x_shift14 = x <<< 14;
    wire signed [40:0] x_shift13 = x <<< 13;
    wire signed [40:0] x_shift12 = x <<< 12;
    wire signed [40:0] x_shift11 = x <<< 11;
    wire signed [40:0] x_shift10 = x <<< 10;
    wire signed [40:0] x_shift9 = x <<< 9;
    wire signed [40:0] x_shift8 = x <<< 8;
    wire signed [40:0] x_shift7 = x <<< 7;
    wire signed [40:0] x_shift4 = x <<< 4;
    wire signed [40:0] x_shift3 = x <<< 3;
    wire signed [40:0] x_shift2 = x <<< 2;
    wire signed [40:0] x_shift0 = x <<< 0;

    // Cross-CSE: shared pattern "+0-"
    wire signed [40:0] _cse_0 = x_shift22 - x_shift20;

    // h0: 000000+0-0+000+00000000000
    assign h0 = (_cse_0 >>> 3) + x_shift15 + x_shift11;

    // h1: 00000+0-0-000-000000000000
    assign h1 = (_cse_0 >>> 2) + -x_shift16 - x_shift12;

    // h2: 00000+0+0-000+000000000000
    assign h2 = x_shift20 + (_cse_0 >>> 4) + x_shift12;

    // h3: 0000+00-0+0+00000000000000
    assign h3 = x_shift21 - x_shift18 + x_shift16 + x_shift14;

    // h4: 0000+0+00+0-00000000000000
    assign h4 = x_shift21 + x_shift19 + (_cse_0 >>> 6);

    // h5: 000+0-0+000+00000000000000
    assign h5 = _cse_0 + x_shift18 + x_shift14;

    // h6: 000+0000-000+00+0000000000
    assign h6 = x_shift22 - x_shift17 + x_shift13 + x_shift10;

    // h7: 000+000+0+00+0000000000000
    assign h7 = x_shift22 + x_shift18 + x_shift16 + x_shift13;

    // h8: 000+00+000000-00+000000000
    assign h8 = x_shift22 + x_shift19 - x_shift12 + x_shift9;

    // h9: 000+00+0-0-000000000000000
    assign h9 = x_shift22 + (_cse_0 >>> 3) + -x_shift15;

    // h10: 000+0000-0-0+0000000000000
    assign h10 = x_shift22 - x_shift17 - x_shift15 + x_shift13;

    // h11: 000+0-000+00+0000000000000
    assign h11 = _cse_0 + x_shift16 + x_shift13;

    // h12: 0000+000+0-0-0000000000000
    assign h12 = x_shift21 + (_cse_0 >>> 5) + -x_shift13;

    // h13: 00000+00000+0+0-0000000000
    assign h13 = x_shift20 + x_shift14 + (_cse_0 >>> 10);

    // h14: 00000000000+0-0000000+000-
    assign h14 = (_cse_0 >>> 8) + x_shift4 - x_shift0;

    // h15: 00000-00+0+0000-0000000000
    assign h15 = -x_shift20 + x_shift17 + x_shift15 - x_shift10;

    // h16: 0000-0+00+00000-0000000000
    assign h16 = -x_shift21 + x_shift19 + x_shift16 - x_shift10;

    // h17: 0000-00+00+000+00000000000
    assign h17 = -x_shift21 + x_shift18 + x_shift15 + x_shift11;

    // h18: 0000-00+0+000000+000000000
    assign h18 = -x_shift21 + x_shift18 + x_shift16 + x_shift9;

    // h19: 0000-0+0+0-000000000000000
    assign h19 = -x_shift21 + x_shift19 + (_cse_0 >>> 5);

    // h20: 00000-000+0-0+000000000000
    assign h20 = -x_shift20 + (_cse_0 >>> 6) + x_shift12;

    // h21: 000000-0+0-0-0000000000000
    assign h21 = -x_shift19 + (_cse_0 >>> 5) + -x_shift13;

    // h22: 00000000+00-00-0+000000000
    assign h22 = x_shift17 - x_shift14 - x_shift11 + x_shift9;

    // h23: 000000+00+0-0000+000000000
    assign h23 = x_shift19 + (_cse_0 >>> 6) + x_shift9;

    // h24: 00000+0-0+0+00000000000000
    assign h24 = (_cse_0 >>> 2) + x_shift16 + x_shift14;

    // h25: 00000+0000-0-0000-00000000
    assign h25 = x_shift20 - x_shift15 - x_shift13 - x_shift8;

    // h26: 00000+000-000-00000000-000
    assign h26 = x_shift20 - x_shift16 - x_shift12 - x_shift3;

    // h27: 00000+0-0+000+000000000000
    assign h27 = (_cse_0 >>> 2) + x_shift16 + x_shift12;

    // h28: 000000+0+000+0+00000000000
    assign h28 = x_shift19 + x_shift17 + x_shift13 + x_shift11;

    // h29: 000000+00-0000+00-00000000
    assign h29 = x_shift19 - x_shift16 + x_shift11 - x_shift8;

    // h30: 0000000+000000+000-0000+00
    assign h30 = x_shift18 + x_shift11 - x_shift7 + x_shift2;

    // h31: 00000000+000-00+0+00000000
    assign h31 = x_shift17 - x_shift13 + x_shift10 + x_shift8;
endmodule
