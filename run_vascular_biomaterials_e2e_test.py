from dotenv import load_dotenv

from biomni.agent import A1

load_dotenv()

agent = A1(
    path='./data',
    llm='claude-sonnet-4-5-20250929',
)

prompt = """Using the biomni.tool.vascular_biomaterials tools, do the following for a
biodegradable small-caliber vascular graft design:

1. Model the degradation of a 3-layer scaffold with target half-lives of 7, 14, and 28
   weeks for the inner, middle, and outer layers respectively (thicknesses 60, 100, and
   150 micrometers), and generate a multilayer design report.
2. Assess the diameter-dependent thrombosis risk for a 3mm diameter graft at 5 mL/s
   blood flow, and compare it against an 8mm diameter graft at the same flow rate.

Save all output files to ./test_output_e2e/
"""

agent.go(prompt)

agent.save_conversation_history("vascular_biomaterials_e2e_result.pdf")
