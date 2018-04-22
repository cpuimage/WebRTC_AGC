/*
 *  Copyright (c) 2011 The WebRTC project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef MODULES_AUDIO_PROCESSING_AGC_LEGACY_ANALOG_AGC_H_
#define MODULES_AUDIO_PROCESSING_AGC_LEGACY_ANALOG_AGC_H_

//#define MIC_LEVEL_FEEDBACK
#ifdef WEBRTC_AGC_DEBUG_DUMP
#include <stdio.h>
#endif


#include <stdint.h>  // NOLINT(build/include)
#include <string.h>

#ifdef WEBRTC_AGC_DEBUG_DUMP
#include <stdio.h>
#endif

#include <stdint.h>  // NOLINT(build/include)

#include <assert.h>
// If you for some reson need to know if DCHECKs are on, test the value of
// RTC_DCHECK_IS_ON. (Test its value, not if it's defined; it'll always be
// defined, to either a true or a false value.)
#if !defined(NDEBUG) || defined(DCHECK_ALWAYS_ON)
#define RTC_DCHECK_IS_ON 1
#else
#define RTC_DCHECK_IS_ON 0
#endif


// C version. Lacks many features compared to the C++ version, but usage
// guidelines are the same.

#define RTC_DCHECK(condition)   assert(condition)

#define RTC_DCHECK_LT(a, b) RTC_DCHECK((a) < (b))
#define RTC_DCHECK_GE(a, b) RTC_DCHECK((a) >= (b))
#define RTC_DCHECK_GT(a, b) RTC_DCHECK((a) > (b))

// the 32 most significant bits of A(19) * B(26) >> 13
#define AGC_MUL32(A, B) (((B) >> 13) * (A) + (((0x00001FFF & (B)) * (A)) >> 13))
// C + the 32 most significant bits of A * B
#define AGC_SCALEDIFF32(A, B, C) \
  ((C) + ((B) >> 16) * (A) + (((0x0000FFFF & (B)) * (A)) >> 16))

// C + the 32 most significant bits of A * B
#define WEBRTC_SPL_SCALEDIFF32(A, B, C) \
    ((C) + ((B) >> 16) * (A) + (((uint32_t)((B) & 0x0000FFFF) * (A)) >> 16))
// allpass filter coefficients.
static const uint16_t kResampleAllpass1[3] = {3284, 24441, 49528};
static const uint16_t kResampleAllpass2[3] = {12199, 37471, 60255};

// Multiply a 32-bit value with a 16-bit value and accumulate to another input:
#define MUL_ACCUM_1(a, b, c) WEBRTC_SPL_SCALEDIFF32(a, b, c)
#define MUL_ACCUM_2(a, b, c) WEBRTC_SPL_SCALEDIFF32(a, b, c)

// Shifting with negative numbers not allowed
// We cannot do casting here due to signed/unsigned problem
#define WEBRTC_SPL_LSHIFT_W32(x, c)     ((x) << (c))
#define WEBRTC_SPL_MIN(A, B)        ((A) < (B) ? (A) : (B))  // Get min value
#define WEBRTC_SPL_MAX(A, B)        ((A) > (B) ? (A) : (B))  // Get max value
#define WEBRTC_SPL_UMUL(a, b) \
    ((uint32_t) ((uint32_t)(a) * (uint32_t)(b)))

#define WEBRTC_SPL_WORD32_MAX       (int32_t)0x7fffffff
#define WEBRTC_SPL_WORD32_MIN       (int32_t)0x80000000
// Shifting with negative numbers allowed
// Positive means left shift
#define WEBRTC_SPL_SHIFT_W32(x, c) ((c) >= 0 ? (x) * (1 << (c)) : (x) >> -(c))

#define WEBRTC_SPL_UMUL_32_16(a, b) \
    ((uint32_t) ((uint32_t)(a) * (uint16_t)(b)))
#define WEBRTC_SPL_MUL_16_U16(a, b) \
    ((int32_t)(int16_t)(a) * (uint16_t)(b))
#define WEBRTC_SPL_ABS_W32(a) \
    (((int32_t)(a) >= 0) ? ((int32_t)(a)) : -((int32_t)(a)))

static __inline int16_t WebRtcSpl_DivW32W16ResW16(int32_t num, int16_t den) {
    // Guard against division with 0
    if (den != 0) {
        return (int16_t) (num / den);
    } else {
        return (int16_t) 0x7FFF;
    }
}

static __inline int32_t WebRtcSpl_DivW32W16(int32_t num, int16_t den) {
    // Guard against division with 0
    if (den != 0) {
        return (int32_t) (num / den);
    } else {
        return (int32_t) 0x7FFFFFFF;
    }
}

// Table used by WebRtcSpl_CountLeadingZeros32_NotBuiltin. For each uint32_t n
// that's a sequence of 0 bits followed by a sequence of 1 bits, the entry at
// index (n * 0x8c0b2891) >> 26 in this table gives the number of zero bits in
// n.
static const int8_t kWebRtcSpl_CountLeadingZeros32_Table[64] = {
        32, 8, 17, -1, -1, 14, -1, -1, -1, 20, -1, -1, -1, 28, -1, 18,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 26, 25, 24,
        4, 11, 23, 31, 3, 7, 10, 16, 22, 30, -1, -1, 2, 6, 13, 9,
        -1, 15, -1, 21, -1, 29, 19, -1, -1, -1, -1, -1, 1, 27, 5, 12,
};

// Don't call this directly except in tests!
static __inline int WebRtcSpl_CountLeadingZeros32_NotBuiltin(uint32_t n) {
    // Normalize n by rounding up to the nearest number that is a sequence of 0
    // bits followed by a sequence of 1 bits. This number has the same number of
    // leading zeros as the original n. There are exactly 33 such values.
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;

    // Multiply the modified n with a constant selected (by exhaustive search)
    // such that each of the 33 possible values of n give a product whose 6 most
    // significant bits are unique. Then look up the answer in the table.
    return kWebRtcSpl_CountLeadingZeros32_Table[(n * 0x8c0b2891) >> 26];
}

// Returns the number of leading zero bits in the argument.
static __inline int WebRtcSpl_CountLeadingZeros32(uint32_t n) {
#ifdef __GNUC__
    return n == 0 ? 32 : __builtin_clz(n);
#else
    return WebRtcSpl_CountLeadingZeros32_NotBuiltin(n);
#endif
}

// Return the number of steps a can be left-shifted without overflow,
// or 0 if a == 0.
static __inline int16_t WebRtcSpl_NormU32(uint32_t a) {
    return a == 0 ? 0 : WebRtcSpl_CountLeadingZeros32(a);
}


static __inline int16_t WebRtcSpl_SatW32ToW16(int32_t value32) {
    int16_t out16 = (int16_t) value32;

    if (value32 > 32767)
        out16 = 32767;
    else if (value32 < -32768)
        out16 = -32768;

    return out16;
}


// Return the number of steps a can be left-shifted without overflow,
// or 0 if a == 0.
static __inline int16_t WebRtcSpl_NormW32(int32_t a) {
    return a == 0 ? 0 : WebRtcSpl_CountLeadingZeros32(a < 0 ? ~a : a) - 1;
}

static __inline int16_t WebRtcSpl_AddSatW16(int16_t a, int16_t b) {
    return WebRtcSpl_SatW32ToW16((int32_t) a + (int32_t) b);
}

static __inline void WebRtcSpl_DownsampleBy2(const int16_t *in, size_t len,
                                             int16_t *out, int32_t *filtState) {
    int32_t tmp1, tmp2, diff, in32, out32;
    size_t i;

    register int32_t state0 = filtState[0];
    register int32_t state1 = filtState[1];
    register int32_t state2 = filtState[2];
    register int32_t state3 = filtState[3];
    register int32_t state4 = filtState[4];
    register int32_t state5 = filtState[5];
    register int32_t state6 = filtState[6];
    register int32_t state7 = filtState[7];

    for (i = (len >> 1); i > 0; i--) {
        // lower allpass filter
        in32 = (int32_t) (*in++) * (1 << 10);
        diff = in32 - state1;
        tmp1 = MUL_ACCUM_1(kResampleAllpass2[0], diff, state0);
        state0 = in32;
        diff = tmp1 - state2;
        tmp2 = MUL_ACCUM_2(kResampleAllpass2[1], diff, state1);
        state1 = tmp1;
        diff = tmp2 - state3;
        state3 = MUL_ACCUM_2(kResampleAllpass2[2], diff, state2);
        state2 = tmp2;

        // upper allpass filter
        in32 = (int32_t) (*in++) * (1 << 10);
        diff = in32 - state5;
        tmp1 = MUL_ACCUM_1(kResampleAllpass1[0], diff, state4);
        state4 = in32;
        diff = tmp1 - state6;
        tmp2 = MUL_ACCUM_1(kResampleAllpass1[1], diff, state5);
        state5 = tmp1;
        diff = tmp2 - state7;
        state7 = MUL_ACCUM_2(kResampleAllpass1[2], diff, state6);
        state6 = tmp2;

        // add two allpass outputs, divide by two and round
        out32 = (state3 + state7 + 1024) >> 11;

        // limit amplitude to prevent wrap-around, and write to output array
        *out++ = WebRtcSpl_SatW32ToW16(out32);
    }

    filtState[0] = state0;
    filtState[1] = state1;
    filtState[2] = state2;
    filtState[3] = state3;
    filtState[4] = state4;
    filtState[5] = state5;
    filtState[6] = state6;
    filtState[7] = state7;
}

typedef struct {
    int32_t downState[8];
    int16_t HPstate;
    int16_t counter;
    int16_t logRatio;           // log( P(active) / P(inactive) ) (Q10)
    int16_t meanLongTerm;       // Q10
    int32_t varianceLongTerm;   // Q8
    int16_t stdLongTerm;        // Q10
    int16_t meanShortTerm;      // Q10
    int32_t varianceShortTerm;  // Q8
    int16_t stdShortTerm;       // Q10
} AgcVad;                     // total = 54 bytes

typedef struct {
    int32_t capacitorSlow;
    int32_t capacitorFast;
    int32_t gain;
    int32_t gainTable[32];
    int16_t gatePrevious;
    int16_t agcMode;
    AgcVad vadNearend;
    AgcVad vadFarend;
#ifdef WEBRTC_AGC_DEBUG_DUMP
    FILE* logFile;
    int frameCounter;
#endif
} DigitalAgc;

int32_t WebRtcAgc_InitDigital(DigitalAgc *digitalAgcInst, int16_t agcMode);

int32_t WebRtcAgc_ProcessDigital(DigitalAgc *digitalAgcInst,
                                 const int16_t *const *inNear,
                                 size_t num_bands,
                                 int16_t *const *out,
                                 uint32_t FS,
                                 int16_t lowLevelSignal);

int32_t WebRtcAgc_AddFarendToDigital(DigitalAgc *digitalAgcInst,
                                     const int16_t *inFar,
                                     size_t nrSamples);

void WebRtcAgc_InitVad(AgcVad *vadInst);

int16_t WebRtcAgc_ProcessVad(AgcVad *vadInst,    // (i) VAD state
                             const int16_t *in,  // (i) Speech signal
                             size_t nrSamples);  // (i) number of samples

int32_t WebRtcAgc_CalculateGainTable(int32_t *gainTable,         // Q16
                                     int16_t compressionGaindB,  // Q0 (in dB)
                                     int16_t targetLevelDbfs,    // Q0 (in dB)
                                     uint8_t limiterEnable,
                                     int16_t analogTarget);

// Errors
#define AGC_UNSPECIFIED_ERROR 18000
#define AGC_UNSUPPORTED_FUNCTION_ERROR 18001
#define AGC_UNINITIALIZED_ERROR 18002
#define AGC_NULL_POINTER_ERROR 18003
#define AGC_BAD_PARAMETER_ERROR 18004

// Warnings
#define AGC_BAD_PARAMETER_WARNING 18050

enum {
    kAgcModeUnchanged,
    kAgcModeAdaptiveAnalog,
    kAgcModeAdaptiveDigital,
    kAgcModeFixedDigital
};

enum {
    kAgcFalse = 0, kAgcTrue
};

typedef struct {
    int16_t targetLevelDbfs;    // default 3 (-3 dBOv)
    int16_t compressionGaindB;  // default 9 dB
    uint8_t limiterEnable;      // default kAgcTrue (on)
} WebRtcAgcConfig;

#if defined(__cplusplus)
extern "C" {
#endif

/*
 * This function analyses the number of samples passed to
 * farend and produces any error code that could arise.
 *
 * Input:
 *      - agcInst           : AGC instance.
 *      - samples           : Number of samples in input vector.
 *
 * Return value:
 *                          :  0 - Normal operation.
 *                          : -1 - Error.
 */
int WebRtcAgc_GetAddFarendError(void *state, size_t samples);

/*
 * This function processes a 10 ms frame of far-end speech to determine
 * if there is active speech. The length of the input speech vector must be
 * given in samples (80 when FS=8000, and 160 when FS=16000, FS=32000 or
 * FS=48000).
 *
 * Input:
 *      - agcInst           : AGC instance.
 *      - inFar             : Far-end input speech vector
 *      - samples           : Number of samples in input vector
 *
 * Return value:
 *                          :  0 - Normal operation.
 *                          : -1 - Error
 */
int WebRtcAgc_AddFarend(void *agcInst, const int16_t *inFar, size_t samples);

/*
 * This function processes a 10 ms frame of microphone speech to determine
 * if there is active speech. The length of the input speech vector must be
 * given in samples (80 when FS=8000, and 160 when FS=16000, FS=32000 or
 * FS=48000). For very low input levels, the input signal is increased in level
 * by multiplying and overwriting the samples in inMic[].
 *
 * This function should be called before any further processing of the
 * near-end microphone signal.
 *
 * Input:
 *      - agcInst           : AGC instance.
 *      - inMic             : Microphone input speech vector for each band
 *      - num_bands         : Number of bands in input vector
 *      - samples           : Number of samples in input vector
 *
 * Return value:
 *                          :  0 - Normal operation.
 *                          : -1 - Error
 */
int WebRtcAgc_AddMic(void *agcInst,
                     int16_t *const *inMic,
                     size_t num_bands,
                     size_t samples);

/*
 * This function replaces the analog microphone with a virtual one.
 * It is a digital gain applied to the input signal and is used in the
 * agcAdaptiveDigital mode where no microphone level is adjustable. The length
 * of the input speech vector must be given in samples (80 when FS=8000, and 160
 * when FS=16000, FS=32000 or FS=48000).
 *
 * Input:
 *      - agcInst           : AGC instance.
 *      - inMic             : Microphone input speech vector for each band
 *      - num_bands         : Number of bands in input vector
 *      - samples           : Number of samples in input vector
 *      - micLevelIn        : Input level of microphone (static)
 *
 * Output:
 *      - inMic             : Microphone output after processing (L band)
 *      - inMic_H           : Microphone output after processing (H band)
 *      - micLevelOut       : Adjusted microphone level after processing
 *
 * Return value:
 *                          :  0 - Normal operation.
 *                          : -1 - Error
 */
int WebRtcAgc_VirtualMic(void *agcInst,
                         int16_t *const *inMic,
                         size_t num_bands,
                         size_t samples,
                         int32_t micLevelIn,
                         int32_t *micLevelOut);

/*
 * This function processes a 10 ms frame and adjusts (normalizes) the gain both
 * analog and digitally. The gain adjustments are done only during active
 * periods of speech. The length of the speech vectors must be given in samples
 * (80 when FS=8000, and 160 when FS=16000, FS=32000 or FS=48000). The echo
 * parameter can be used to ensure the AGC will not adjust upward in the
 * presence of echo.
 *
 * This function should be called after processing the near-end microphone
 * signal, in any case after any echo cancellation.
 *
 * Input:
 *      - agcInst           : AGC instance
 *      - inNear            : Near-end input speech vector for each band
 *      - num_bands         : Number of bands in input/output vector
 *      - samples           : Number of samples in input/output vector
 *      - inMicLevel        : Current microphone volume level
 *      - echo              : Set to 0 if the signal passed to add_mic is
 *                            almost certainly free of echo; otherwise set
 *                            to 1. If you have no information regarding echo
 *                            set to 0.
 *
 * Output:
 *      - outMicLevel       : Adjusted microphone volume level
 *      - out               : Gain-adjusted near-end speech vector
 *                          : May be the same vector as the input.
 *      - saturationWarning : A returned value of 1 indicates a saturation event
 *                            has occurred and the volume cannot be further
 *                            reduced. Otherwise will be set to 0.
 *
 * Return value:
 *                          :  0 - Normal operation.
 *                          : -1 - Error
 */
int WebRtcAgc_Process(void *agcInst,
                      const int16_t *const *inNear,
                      size_t num_bands,
                      size_t samples,
                      int16_t *const *out,
                      int32_t inMicLevel,
                      int32_t *outMicLevel,
                      int16_t echo,
                      uint8_t *saturationWarning);

/*
 * This function sets the config parameters (targetLevelDbfs,
 * compressionGaindB and limiterEnable).
 *
 * Input:
 *      - agcInst           : AGC instance
 *      - config            : config struct
 *
 * Output:
 *
 * Return value:
 *                          :  0 - Normal operation.
 *                          : -1 - Error
 */
int WebRtcAgc_set_config(void *agcInst, WebRtcAgcConfig config);

/*
 * This function returns the config parameters (targetLevelDbfs,
 * compressionGaindB and limiterEnable).
 *
 * Input:
 *      - agcInst           : AGC instance
 *
 * Output:
 *      - config            : config struct
 *
 * Return value:
 *                          :  0 - Normal operation.
 *                          : -1 - Error
 */
int WebRtcAgc_get_config(void *agcInst, WebRtcAgcConfig *config);

/*
 * This function creates and returns an AGC instance, which will contain the
 * state information for one (duplex) channel.
 */
void *WebRtcAgc_Create(void);

/*
 * This function frees the AGC instance created at the beginning.
 *
 * Input:
 *      - agcInst           : AGC instance.
 */
void WebRtcAgc_Free(void *agcInst);

/*
 * This function initializes an AGC instance.
 *
 * Input:
 *      - agcInst           : AGC instance.
 *      - minLevel          : Minimum possible mic level
 *      - maxLevel          : Maximum possible mic level
 *      - agcMode           : 0 - Unchanged
 *                          : 1 - Adaptive Analog Automatic Gain Control -3dBOv
 *                          : 2 - Adaptive Digital Automatic Gain Control -3dBOv
 *                          : 3 - Fixed Digital Gain 0dB
 *      - fs                : Sampling frequency
 *
 * Return value             :  0 - Ok
 *                            -1 - Error
 */
int WebRtcAgc_Init(void *agcInst,
                   int32_t minLevel,
                   int32_t maxLevel,
                   int16_t agcMode,
                   uint32_t fs);

#if defined(__cplusplus)
}
#endif

/* Analog Automatic Gain Control variables:
 * Constant declarations (inner limits inside which no changes are done)
 * In the beginning the range is narrower to widen as soon as the measure
 * 'Rxx160_LP' is inside it. Currently the starting limits are -22.2+/-1dBm0
 * and the final limits -22.2+/-2.5dBm0. These levels makes the speech signal
 * go towards -25.4dBm0 (-31.4dBov). Tuned with wbfile-31.4dBov.pcm
 * The limits are created by running the AGC with a file having the desired
 * signal level and thereafter plotting Rxx160_LP in the dBm0-domain defined
 * by out=10*log10(in/260537279.7); Set the target level to the average level
 * of our measure Rxx160_LP. Remember that the levels are in blocks of 16 in
 * Q(-7). (Example matlab code: round(db2pow(-21.2)*16/2^7) )
 */
#define RXX_BUFFER_LEN 10

static const int16_t kMsecSpeechInner = 520;
static const int16_t kMsecSpeechOuter = 340;

static const int16_t kNormalVadThreshold = 400;

static const int16_t kAlphaShortTerm = 6;  // 1 >> 6 = 0.0156
static const int16_t kAlphaLongTerm = 10;  // 1 >> 10 = 0.000977

typedef struct {
    // Configurable parameters/variables
    uint32_t fs;                // Sampling frequency
    int16_t compressionGaindB;  // Fixed gain level in dB
    int16_t targetLevelDbfs;    // Target level in -dBfs of envelope (default -3)
    int16_t agcMode;            // Hard coded mode (adaptAna/adaptDig/fixedDig)
    uint8_t limiterEnable;      // Enabling limiter (on/off (default off))
    WebRtcAgcConfig defaultConfig;
    WebRtcAgcConfig usedConfig;

    // General variables
    int16_t initFlag;
    int16_t lastError;

    // Target level parameters
    // Based on the above: analogTargetLevel = round((32767*10^(-22/20))^2*16/2^7)
    int32_t analogTargetLevel;    // = RXX_BUFFER_LEN * 846805;       -22 dBfs
    int32_t startUpperLimit;      // = RXX_BUFFER_LEN * 1066064;      -21 dBfs
    int32_t startLowerLimit;      // = RXX_BUFFER_LEN * 672641;       -23 dBfs
    int32_t upperPrimaryLimit;    // = RXX_BUFFER_LEN * 1342095;      -20 dBfs
    int32_t lowerPrimaryLimit;    // = RXX_BUFFER_LEN * 534298;       -24 dBfs
    int32_t upperSecondaryLimit;  // = RXX_BUFFER_LEN * 2677832;      -17 dBfs
    int32_t lowerSecondaryLimit;  // = RXX_BUFFER_LEN * 267783;       -27 dBfs
    uint16_t targetIdx;           // Table index for corresponding target level
#ifdef MIC_LEVEL_FEEDBACK
    uint16_t targetIdxOffset;  // Table index offset for level compensation
#endif
    int16_t analogTarget;  // Digital reference level in ENV scale

    // Analog AGC specific variables
    int32_t filterState[8];  // For downsampling wb to nb
    int32_t upperLimit;      // Upper limit for mic energy
    int32_t lowerLimit;      // Lower limit for mic energy
    int32_t Rxx160w32;       // Average energy for one frame
    int32_t Rxx16_LPw32;     // Low pass filtered subframe energies
    int32_t Rxx160_LPw32;    // Low pass filtered frame energies
    int32_t Rxx16_LPw32Max;  // Keeps track of largest energy subframe
    int32_t Rxx16_vectorw32[RXX_BUFFER_LEN];  // Array with subframe energies
    int32_t Rxx16w32_array[2][5];  // Energy values of microphone signal
    int32_t env[2][10];            // Envelope values of subframes

    int16_t Rxx16pos;          // Current position in the Rxx16_vectorw32
    int16_t envSum;            // Filtered scaled envelope in subframes
    int16_t vadThreshold;      // Threshold for VAD decision
    int16_t inActive;          // Inactive time in milliseconds
    int16_t msTooLow;          // Milliseconds of speech at a too low level
    int16_t msTooHigh;         // Milliseconds of speech at a too high level
    int16_t changeToSlowMode;  // Change to slow mode after some time at target
    int16_t firstCall;         // First call to the process-function
    int16_t msZero;            // Milliseconds of zero input
    int16_t msecSpeechOuterChange;  // Min ms of speech between volume changes
    int16_t msecSpeechInnerChange;  // Min ms of speech between volume changes
    int16_t activeSpeech;           // Milliseconds of active speech
    int16_t muteGuardMs;            // Counter to prevent mute action
    int16_t inQueue;                // 10 ms batch indicator

    // Microphone level variables
    int32_t micRef;         // Remember ref. mic level for virtual mic
    uint16_t gainTableIdx;  // Current position in virtual gain table
    int32_t micGainIdx;     // Gain index of mic level to increase slowly
    int32_t micVol;         // Remember volume between frames
    int32_t maxLevel;       // Max possible vol level, incl dig gain
    int32_t maxAnalog;      // Maximum possible analog volume level
    int32_t maxInit;        // Initial value of "max"
    int32_t minLevel;       // Minimum possible volume level
    int32_t minOutput;      // Minimum output volume level
    int32_t zeroCtrlMax;    // Remember max gain => don't amp low input
    int32_t lastInMicLevel;

    int16_t scale;  // Scale factor for internal volume levels
#ifdef MIC_LEVEL_FEEDBACK
    int16_t numBlocksMicLvlSat;
    uint8_t micLvlSat;
#endif
    // Structs for VAD and digital_agc
    AgcVad vadMic;
    DigitalAgc digitalAgc;

#ifdef WEBRTC_AGC_DEBUG_DUMP
    FILE* fpt;
    FILE* agcLog;
    int32_t fcount;
#endif

    int16_t lowLevelSignal;
} LegacyAgc;

#endif  // MODULES_AUDIO_PROCESSING_AGC_LEGACY_ANALOG_AGC_H_
