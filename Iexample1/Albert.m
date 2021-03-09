def createDIRK_OrderThreeStagesFour_WSO2(self):
        # See file:///Users/albertcosta/Downloads/2020_Book_SpectralAndHighOrderMethodsFor.pdf pag 448
        rungeKuttaParameters = RK.RungeKuttaParameters()

        def CoefA():
            theCoefA = np.zeros((4, 4), dtype='float')
            theCoefA[0, 0] = 0.01900072890

            theCoefA[1, 0] = 0.40434605601
            theCoefA[1, 1] = 0.38435717512

            theCoefA[2, 0] = 0.06487908412
            theCoefA[2, 1] = -0.16389640295
            theCoefA[2, 2] = 0.51545231222

            theCoefA[3, 0] = 0.02343549374
            theCoefA[3, 1] = -0.41207877888
            theCoefA[3, 2] = 0.96661161281
            theCoefA[3, 3] = 0.42203167233

            nOStages = theCoefA.shape[0]
            stagesArray = np.arange(nOStages)

            theCoefInvA = np.tril(np.linalg.inv(theCoefA))
            theSumCoefInvA = theCoefInvA[stagesArray, :].sum(1)

            return theCoefA, theCoefInvA, theSumCoefInvA

        def CoefB():
            theCoefB = np.array([0.02343549374, -0.41207877888, 0.96661161281, 0.42203167233])
            return theCoefB

        def CoefC():
            theCoefC = np.array([0.01900072890, 0.78870323114, 0.41643499339, 1.0])
            return theCoefC

        rungeKuttaParameters.theCoefA = CoefA
        rungeKuttaParameters.theCoefB = CoefB
        rungeKuttaParameters.theCoefC = CoefC

        return rungeKuttaParameters