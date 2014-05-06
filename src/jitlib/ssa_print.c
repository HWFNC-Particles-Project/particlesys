#include <ssa_def.h>
#include <stdio.h>
#include <inttypes.h>

static void print_flags(uint64_t op) {
    printf("        ");
    if(op & DEAD_BIT) {
        printf(" DEAD");
    }
}

static void print_argument(ssa_block *block, uint64_t arg0) {
    uint64_t arg1 = ssa_get_index(block, arg0);
    if(arg0 == arg1) {
        printf("s%" PRIu64, arg0);
    } else {
        printf("(s%"PRIu64"->s%"PRIu64")", arg0, arg1);
    }
}

void ssa_print_op(ssa_block *block, size_t i) {
   uint64_t op = block->buffer[i];
   switch(GET_OP(op)) {
        case SSA_STORE:
            printf("[%"PRIu64"] = ", GET_ARG2(op));
            print_argument(block, GET_ARG1(op));
            break;
        case SSA_ADD:
            printf("s%zu = ", i);
            print_argument(block, GET_ARG1(op));
            printf(" + ");
            print_argument(block, GET_ARG2(op));
            break;
        case SSA_SUB:
            printf("s%zu = ", i);
            print_argument(block, GET_ARG1(op));
            printf(" - ");
            print_argument(block, GET_ARG2(op));
            break;
        case SSA_MUL:
            printf("s%zu = ", i);
            print_argument(block, GET_ARG1(op));
            printf(" * ");
            print_argument(block, GET_ARG2(op));
            break;
        case SSA_DIV:
            printf("s%zu = ", i);
            print_argument(block, GET_ARG1(op));
            printf(" / ");
            print_argument(block, GET_ARG2(op));
            break;
        case SSA_MIN:
            printf("s%zu = min(", i);
            print_argument(block, GET_ARG1(op));
            printf(", ");
            print_argument(block, GET_ARG2(op));
            printf(")");
            break;
        case SSA_MAX:
            printf("s%zu = max(", i);
            print_argument(block, GET_ARG1(op));
            printf(", ");
            print_argument(block, GET_ARG2(op));
            printf(")");
            break;
        case SSA_AND:
            printf("s%zu = ", i);
            print_argument(block, GET_ARG1(op));
            printf(" & ");
            print_argument(block, GET_ARG2(op));
            break;
        case SSA_ANDN:
            printf("s%zu = ", i);
            print_argument(block, GET_ARG1(op));
            printf(" &! ");
            print_argument(block, GET_ARG2(op));
            break;
        case SSA_OR:
            printf("s%zu = ", i);
            print_argument(block, GET_ARG1(op));
            printf(" | ");
            print_argument(block, GET_ARG2(op));
            break;
        case SSA_XOR:
            printf("s%zu = ", i);
            print_argument(block, GET_ARG1(op));
            printf(" ^ ");
            print_argument(block, GET_ARG2(op));
            break;

        case SSA_CMPEQ:
            printf("s%zu = ", i);
            print_argument(block, GET_ARG1(op));
            printf(" == ");
            print_argument(block, GET_ARG2(op));
            break;
        case SSA_CMPLE:
            printf("s%zu = ", i);
            print_argument(block, GET_ARG1(op));
            printf(" <= ");
            print_argument(block, GET_ARG2(op));
            break;
        case SSA_CMPLT:
            printf("s%zu = ", i);
            print_argument(block, GET_ARG1(op));
            printf(" < ");
            print_argument(block, GET_ARG2(op));
            break;

        case SSA_CMPNEQ:
            printf("s%zu = ", i);
            print_argument(block, GET_ARG1(op));
            printf(" != ");
            print_argument(block, GET_ARG2(op));
            break;
        case SSA_CMPNLE:
            printf("s%zu = ", i);
            print_argument(block, GET_ARG1(op));
            printf(" > ");
            print_argument(block, GET_ARG2(op));
            break;
        case SSA_CMPNLT:
            printf("s%zu = ", i);
            print_argument(block, GET_ARG1(op));
            printf(" >= ");
            print_argument(block, GET_ARG2(op));
            break;

        case SSA_LOAD:
            printf("s%zu = [%"PRIu64"]", i, GET_ARG1(op));
            break;
        case SSA_CONST:
            printf("s%zu = %f", i, GET_CONST(op));
            break;
        case SSA_SQRT:
            printf("s%zu = sqrt(", i);
            print_argument(block, GET_ARG1(op));
            printf(")");
            break;
        case SSA_RSQRT:
            printf("s%zu = rsqrt(", i);
            print_argument(block, GET_ARG1(op));
            printf(")");
            break;
        case SSA_RCP:
            printf("s%zu = rcp(", i);
            print_argument(block, GET_ARG1(op));
            printf(")");
            break;

        case SSA_BLEND:
            printf("s%zu = ", i);
            print_argument(block, GET_ARG1(op));
            printf(" ? ");
            print_argument(block, GET_ARG2(op));
            printf(" : ");
            print_argument(block, GET_ARG3(op));
            break;

        default:
            printf("unknown");
            break;
    }
}

void ssa_print(ssa_block *block) {
    for(size_t i = 0;i<block->size;++i) {
        ssa_print_op(block, i);
        print_flags(block->buffer[i]);
        printf("\n");
    }
}
