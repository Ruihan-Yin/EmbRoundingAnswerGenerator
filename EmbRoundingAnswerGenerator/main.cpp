#include <cfenv>
#include <iostream>
#include <immintrin.h>  // Include Intel Intrinsics header
#include <type_traits>
#include <fstream>
#include <string>
#include <iomanip>


__m512d PerformOperationUnary(double val, const int& roundingMode, const std::string& oper)
{
    __m512d valVec = _mm512_set1_pd(val);

    if (oper == std::string("Sqrt"))
    {
        if (roundingMode == 0x09)
        {
            return _mm512_sqrt_round_pd(valVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm512_sqrt_round_pd(valVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm512_sqrt_round_pd(valVec, 0x0B);
        }
    }
    else
    {
        // Handle other operations as needed
        std::cerr << "Unsupported operation: " << oper << std::endl;
        return _mm512_setzero_pd();  // Return default value for unsupported operations
    }
}

__m512 PerformOperationUnary(float val, const int& roundingMode, const std::string& oper)
{
    __m512 valVec = _mm512_set1_ps(val);

    if (oper == std::string("Sqrt"))
    {
        if (roundingMode == 0x09)
        {
            return _mm512_sqrt_round_ps(valVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm512_sqrt_round_ps(valVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm512_sqrt_round_ps(valVec, 0x0B);
        }
    }
    else
    {
        // Handle other operations as needed
        std::cerr << "Unsupported operation: " << oper << std::endl;
        return _mm512_setzero_ps();  // Return default value for unsupported operations
    }
}

__m512d PerformOperationBinary(double left, double right, const int& roundingMode, const std::string& oper)
{
    __m512d leftVec = _mm512_set1_pd(left);
    __m512d rightVec = _mm512_set1_pd(right);
    if (oper == std::string("Add"))
    {
        if (roundingMode == 0x09)
        {
            return _mm512_add_round_pd(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm512_add_round_pd(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm512_add_round_pd(leftVec, rightVec, 0x0B);
        }
    }
    else if (oper == std::string("Subtract"))
    {
        if (roundingMode == 0x09)
        {
            return _mm512_sub_round_pd(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm512_sub_round_pd(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm512_sub_round_pd(leftVec, rightVec, 0x0B);
        }
    }
    else if (oper == std::string("Divide"))
    {
        if (roundingMode == 0x09)
        {
            return _mm512_div_round_pd(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm512_div_round_pd(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm512_div_round_pd(leftVec, rightVec, 0x0B);
        }
    }
    else if (oper == std::string("Multiply"))
    {
        if (roundingMode == 0x09)
        {
            return _mm512_mul_round_pd(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm512_mul_round_pd(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm512_mul_round_pd(leftVec, rightVec, 0x0B);
        }
    }
    else if (oper == std::string("Scale"))
    {
        if (roundingMode == 0x09)
        {
            return _mm512_scalef_round_pd(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm512_scalef_round_pd(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm512_scalef_round_pd(leftVec, rightVec, 0x0B);
        }
    }
    else
    {
        // Handle other operations as needed
        std::cerr << "Unsupported operation: " << oper << std::endl;
        return _mm512_setzero_pd();  // Return default value for unsupported operations
    }
}

__m512 PerformOperationBinary(float left, float right, const int& roundingMode, const std::string& oper)
{

    __m512 leftVec = _mm512_set1_ps(left);
    __m512 rightVec = _mm512_set1_ps(right);
    if (oper == std::string("Add"))
    {

        if (roundingMode == 0x09)
        {
            return _mm512_add_round_ps(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm512_add_round_ps(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm512_add_round_ps(leftVec, rightVec, 0x0B);
        }
    }
    else if (oper == std::string("Subtract"))
    {
        if (roundingMode == 0x09)
        {
            return _mm512_sub_round_ps(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm512_sub_round_ps(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm512_sub_round_ps(leftVec, rightVec, 0x0B);
        }
    }
    else if (oper == std::string("Divide"))
    {
        if (roundingMode == 0x09)
        {
            return _mm512_div_round_ps(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm512_div_round_ps(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm512_div_round_ps(leftVec, rightVec, 0x0B);
        }
    }
    else if (oper == std::string("Multiply"))
    {
        if (roundingMode == 0x09)
        {
            return _mm512_mul_round_ps(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm512_mul_round_ps(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm512_mul_round_ps(leftVec, rightVec, 0x0B);
        }
    }
    else if (oper == std::string("Scale"))
    {
        if (roundingMode == 0x09)
        {
            return _mm512_scalef_round_ps(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm512_scalef_round_ps(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm512_scalef_round_ps(leftVec, rightVec, 0x0B);
        }
    }
    else
    {
        // Handle other operations as needed
        std::cerr << "Unsupported operation: " << oper << std::endl;
        return _mm512_setzero_ps();  // Return default value for unsupported operations
    }
}

__m128d PerformOperationScalarBinary(double left, double right, const int& roundingMode, const std::string& oper)
{

    __m128d leftVec = _mm_set1_pd(left);
    __m128d rightVec = _mm_set1_pd(right);
    if (oper == std::string("AddScalar"))
    {
        if (roundingMode == 0x09)
        {
            return _mm_add_round_sd(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm_add_round_sd(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm_add_round_sd(leftVec, rightVec, 0x0B);
        }
    }
    else if (oper == std::string("SubtractScalar"))
    {
        if (roundingMode == 0x09)
        {
            return _mm_sub_round_sd(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm_sub_round_sd(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm_sub_round_sd(leftVec, rightVec, 0x0B);
        }
    }
    else if (oper == std::string("DivideScalar"))
    {
        if (roundingMode == 0x09)
        {
            return _mm_div_round_sd(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm_div_round_sd(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm_div_round_sd(leftVec, rightVec, 0x0B);
        }
    }
    else if (oper == std::string("MultiplyScalar"))
    {
        if (roundingMode == 0x09)
        {
            return _mm_mul_round_sd(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm_mul_round_sd(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm_mul_round_sd(leftVec, rightVec, 0x0B);
        }
    }
    else if (oper == std::string("SqrtScalar"))
    {
        if (roundingMode == 0x09)
        {
            return _mm_sqrt_round_sd(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm_sqrt_round_sd(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm_sqrt_round_sd(leftVec, rightVec, 0x0B);
        }
    }
    else
    {
        // Handle other operations as needed
        std::cerr << "Unsupported operation: " << oper << std::endl;
        return _mm_setzero_pd();  // Return default value for unsupported operations
    }
}

__m128 PerformOperationScalarBinary(float left, float right, const int& roundingMode, const std::string& oper)
{
    __m128 leftVec = _mm_set1_ps(left);
    __m128 rightVec = _mm_set1_ps(right);
    if (oper == std::string("AddScalar"))
    {
        if (roundingMode == 0x09)
        {
            return _mm_add_round_ss(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm_add_round_ss(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm_add_round_ss(leftVec, rightVec, 0x0B);
        }
    }
    else if (oper == std::string("SubtractScalar"))
    {
        if (roundingMode == 0x09)
        {
            return _mm_sub_round_ss(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm_sub_round_ss(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm_sub_round_ss(leftVec, rightVec, 0x0B);
        }
    }
    else if (oper == std::string("DivideScalar"))
    {
        if (roundingMode == 0x09)
        {
            return _mm_div_round_ss(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm_div_round_ss(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm_div_round_ss(leftVec, rightVec, 0x0B);
        }
    }
    else if (oper == std::string("MultiplyScalar"))
    {
        if (roundingMode == 0x09)
        {
            return _mm_mul_round_ss(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm_mul_round_ss(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm_mul_round_ss(leftVec, rightVec, 0x0B);
        }
    }
    else if (oper == std::string("SqrtScalar"))
    {
        if (roundingMode == 0x09)
        {
            return _mm_sqrt_round_ss(leftVec, rightVec, 0x09);
        }
        else if (roundingMode == 0x0A)
        {
            return _mm_sqrt_round_ss(leftVec, rightVec, 0x0A);
        }
        else if (roundingMode == 0x0B)
        {
            return _mm_sqrt_round_ss(leftVec, rightVec, 0x0B);
        }
    }
    else
    {
        // Handle other operations as needed
        std::cerr << "Unsupported operation: " << oper << std::endl;
        return _mm_setzero_ps();  // Return default value for unsupported operations
    }
}

std::string parseRoundingMode(int roundingMode)
{
    switch (roundingMode)
    {
    case _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC:
        return "ToNegativeInfinity";

    case _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC:
        return "ToPositiveInfinity";

    case _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC:
        return "ToZero";
    default:
        return "#error:wrong rounding mode";
    }
}

int main() {
    std::ofstream outputFile("output.txt");

    if (!outputFile.is_open()) {
        std::cerr << "Unable to open output file." << std::endl;
        return 1;
    }

    unsigned int floatOutputArr[16];
    float floatOutputArrLen8[8];
    double left = 0.05, right = 0.45, val = 29.37;
    const int roundingModes[] = { _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC , _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC , _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC };
    const std::string opers[] = { "Add", "Divide", "Multiply", "Subtract", "Scale" };
    const std::string scalarOpers[] = { "AddScalar","DivideScalar", "MultiplyScalar", "SubtractScalar", "SqrtScalar" };
    const std::string unaryOpers[] = { "Sqrt" };

    for (int i = 0; i < sizeof(unaryOpers) / sizeof(unaryOpers[0]); i++)
    {
        for (int j = 0; j < sizeof(roundingModes) / sizeof(roundingModes[0]); j++)
        {
            __m512d doubleRes = PerformOperationUnary(val, roundingModes[j], unaryOpers[i]);
            __m512  floatRes = PerformOperationUnary((float)val, roundingModes[j], unaryOpers[i]);

            __m512i doubleHex = _mm512_castpd_si512(doubleRes);
            __m512i floatHex = _mm512_castps_si512(floatRes);

            alignas(64) uint64_t doubleHexArr[8];
            alignas(32) uint32_t floatHexArr[16];
            _mm512_storeu_si512((__m512i*)doubleHexArr, doubleHex);
            _mm512_storeu_si512((__m512i*)floatHexArr, floatHex);
            outputFile << "{(\"" << "Double" << "\", " << "\"" << unaryOpers[i] << "\", \"" << parseRoundingMode(roundingModes[j]) << "\")" << ", new ulong[] {";
            for (int d = 0; d < 8; d++)
            {
                outputFile << "0x" << std::hex << doubleHexArr[d];
                if (d != 7) outputFile << ", ";
            }
            outputFile << "}}," << std::endl;

            outputFile << "{(\"" << "Single" << "\", " << "\"" << unaryOpers[i] << "\", \"" << parseRoundingMode(roundingModes[j]) << "\")" << ", new ulong[] {";
            for (int f = 0; f < 16; f++)
            {
                outputFile << "0x" << std::hex << floatHexArr[f];
                if (f != 15) outputFile << ", ";
            }
            outputFile << "}}," << std::endl;
        }
    }

    for (int i = 0; i < sizeof(scalarOpers) / sizeof(scalarOpers[0]); ++i)
    {
        for (int j = 0; j < sizeof(roundingModes) / sizeof(roundingModes[0]); j++)
        {
            __m128d doubleScalarRes = PerformOperationScalarBinary(left, right, roundingModes[j], scalarOpers[i]);
            __m128  floatScalarRes = PerformOperationScalarBinary((float)left, (float)right, roundingModes[j], scalarOpers[i]);

            __m128i doubleSclarHex = _mm_castpd_si128(doubleScalarRes);
            __m128i floatScalarHex = _mm_castps_si128(floatScalarRes);


            alignas(64) uint64_t doubleScalarHexArr[2];
            alignas(32) uint32_t floatScalarHexArr[4];

            _mm_storeu_si128((__m128i*)doubleScalarHexArr, doubleSclarHex);
            _mm_storeu_si128((__m128i*)floatScalarHexArr, floatScalarHex);

            outputFile << "{(\"" << "Double" << "\", " << "\"" << scalarOpers[i] << "\", \"" << parseRoundingMode(roundingModes[j]) << "\")" << ", new ulong[] {";
            for (int s = 0; s < 2; s++)
            {
                outputFile << "0x" << std::hex << doubleScalarHexArr[s];
                if (s != 1) outputFile << ", ";
            }
            outputFile << "}}," << std::endl;

            outputFile << "{(\"" << "Single" << "\", " << "\"" << scalarOpers[i] << "\", \"" << parseRoundingMode(roundingModes[j]) << "\")" << ", new ulong[] {";
            for (int s = 0; s < 4; s++)
            {
                outputFile << "0x" << std::hex << floatScalarHexArr[s];
                if (s != 3) outputFile << ", ";
            }
            outputFile << "}}," << std::endl;
        }
    }

    for (int i = 0; i < sizeof(opers) / sizeof(opers[0]); ++i)
    {
        for (int j = 0; j < sizeof(roundingModes) / sizeof(roundingModes[0]); j++)
        {
            __m512d doubleRes = PerformOperationBinary(left, right, roundingModes[j], opers[i]);
            __m512  floatRes = PerformOperationBinary((float)left, (float)right, roundingModes[j], opers[i]);

            __m512i doubleHex = _mm512_castpd_si512(doubleRes);
            __m512i floatHex = _mm512_castps_si512(floatRes);

            alignas(64) uint64_t doubleHexArr[8];
            alignas(32) uint32_t floatHexArr[16];

            _mm512_storeu_si512((__m512i*)doubleHexArr, doubleHex);
            _mm512_storeu_si512((__m512i*)floatHexArr, floatHex);

            outputFile << "{(\"" << "Double" << "\", " << "\"" << opers[i] << "\", \"" << parseRoundingMode(roundingModes[j]) << "\")" << ", new ulong[] {";
            for (int d = 0; d < 8; d++)
            {
                outputFile << "0x" << std::hex << doubleHexArr[d];
                if (d != 7) outputFile << ", ";
            }
            outputFile << "}}," << std::endl;

            outputFile << "{(\"" << "Single" << "\", " << "\"" << opers[i] << "\", \"" << parseRoundingMode(roundingModes[j]) << "\")" << ", new ulong[] {";
            for (int f = 0; f < 16; f++)
            {
                outputFile << "0x" << std::hex << floatHexArr[f];
                if (f != 15) outputFile << ", ";
            }
            outputFile << "}}," << std::endl;
        }
    }

    // Close the file
    outputFile.close();

    return 0;
}